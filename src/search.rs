// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use genevo::operator::prelude::*;
use genevo::prelude::*;
use genevo::types::fmt::Display;

use gosh::gchemol::prelude::*;
use gosh::gchemol::{Atom, Molecule, io};

use crate::common::*;
use crate::config::Config;
// imports:1 ends here

// base

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*base][base:1]]
type MoleculeGenome = Vec<(usize, f64, f64, f64)>;

/// FIXME: bond connectivity will be lost
fn molecule_from_genome(genome: &MoleculeGenome) -> Molecule {
    let mut mol = Molecule::new("genome");
    for (n, x, y, z) in genome {
        let a = Atom::build().element(*n).position(*x, *y, *z).finish();
        mol.add_atom(a);
    }

    mol.rebond();
    mol
}

fn molecule_to_genome(mol: &Molecule) -> MoleculeGenome {
    let mut g = vec![];

    for a in mol.atoms() {
        let n = a.number();
        let [x, y, z] = a.position();
        g.push((n, x, y, z));
    }

    g
}
// base:1 ends here

// model

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*model][model:1]]
use gosh::models::*;

/// Evaluate energy of molecule in current geometry
fn get_energy(mol: &Molecule, runfile: &str) -> Result<f64> {
    debug!("calculate single point energy using bbm model ...");

    let mut bbm = BlackBox::from_dir(runfile);
    let mr = bbm.compute(mol)?;

    if let Some(energy) = mr.energy {
        debug!("sp energy = {:-12.5}", energy);
        return Ok(energy);
    } else {
        bail!("no energy record found in the output!");
    }
}

/// Return optimized geometry
fn get_optimized_molecule(mol: &Molecule, runfile: &str) -> Result<Molecule> {
    debug!("optimize geometry using bbm model ...");

    // avoid bad geometry
    let mut mol = mol.clone();
    mol.rebond();
    mol.clean()?;

    let mut bbm = BlackBox::from_dir(runfile);

    match bbm.compute(&mol) {
        Ok(mr) => {
            if let Some(energy) = mr.energy {
                info!("opt energy = {:-12.5}", energy);
                if let Some(mol) = mr.molecule {
                    return Ok(mol);
                } else {
                    bail!("no molecule record found in bbm results");
                }
            } else {
                bail!("no energy record found in the output!");
            }
        }

        Err(e) => {
            error!("opt failed, use initial geometry instead.");
            error!("{:?}", e);

            Ok(mol.to_owned())
        }
    }
}
// model:1 ends here

// fitness

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*fitness][fitness:1]]
// energy => fitness
fn calc_fitness(energy: f64, eref: f64) -> u32 {
    let temperature = 50000.;
    let value = (energy - eref) * 96.;
    let fitness = (-1.0 * value / (temperature * 0.0083145)).exp();

    (fitness * 1000.) as u32
}

#[test]
fn test_calc_fit() {
    let x = -811.22008;
    let x = calc_fitness(x, -811.2);
    dbg!(x);
}

#[derive(Clone, Debug, PartialEq)]
struct ClusterFitnessEvaluator<'a> {
    temperature: f64,
    runfile: &'a str,
    eref: f64,
}

impl<'a> ClusterFitnessEvaluator<'a> {
    fn new(runfile: &'a str, eref: f64) -> Self {
        ClusterFitnessEvaluator {
            runfile: &runfile,
            eref: eref,
            temperature: 20000.0,
        }
    }
}

impl<'a> FitnessFunction<MoleculeGenome, u32> for ClusterFitnessEvaluator<'a> {
    fn fitness_of(&self, individual: &MoleculeGenome) -> u32 {
        debug!("evaluate fitness ...");
        let mol = molecule_from_genome(&individual);

        // mapping energy to fitness
        let fitness = {
            if let Ok(energy) = get_energy(&mol, &self.runfile) {
                let fitness = calc_fitness(energy, self.eref);
                println!("final energy = {:-12.5}, fitness= {}", energy, fitness);
                fitness
            } else {
                warn!("no energy for {}", &mol.name);
                self.lowest_possible_fitness()
            }
        };

        fitness
    }

    fn average(&self, fitness_values: &[u32]) -> u32 {
        (fitness_values.iter().sum::<u32>() / fitness_values.len() as u32)
    }

    fn highest_possible_fitness(&self) -> u32 {
        std::u32::MAX
    }

    fn lowest_possible_fitness(&self) -> u32 {
        0
    }
}
// fitness:1 ends here

// operators

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*operators][operators:1]]
use rand::prelude::*;

use genevo::genetic::{Children, Genotype, Parents};
use genevo::operator::prelude::*;
use genevo::operator::{CrossoverOp, GeneticOperator, MutationOp};
use genevo::prelude::*;
use genevo::random::Rng;

use crate::kick;
use spdkit;

/// create a population with n individuals.
/// initial molecule will be read from `molfile`
fn get_population_and_eref(config: &Config) -> (Population<MoleculeGenome>, f64) {
    info!("create initial population ..");
    let mol = Molecule::from_file(&config.molfile).expect("mol2");
    let mut individuals = vec![];
    let n = config.search.population_size;

    // get reference energy
    let mut eref = std::f64::MAX;
    for i in 0..n {
        let mol = kick(&mol).expect("kick mol");
        let mol = get_optimized_molecule(&mol, &config.runfile_opt).expect("opt mol");
        let g = molecule_to_genome(&mol);
        let e = get_energy(&mol, &config.runfile_sp).expect("sp energy");
        if e < eref {
            eref = e;
        }
        individuals.push(g);
    }

    (Population::with_individuals(individuals), eref)
}

#[derive(Clone, Debug, PartialEq)]
pub struct KickMutator<'a> {
    mutation_rate: f64,
    config: &'a Config,
}

impl<'a> KickMutator<'a> {
    pub fn mutate_molecule(&self, mol: &mut Molecule, mutate: bool) -> Result<Molecule> {
        use rand::seq::SliceRandom;

        if mutate {
            let nbonds = mol.nbonds();
            if nbonds > 0 {
                // remove bonds randomly to break molecule into parts
                let mut edges: Vec<_> = mol.bonds().map(|b| b.index()).collect();
                let mut rng = thread_rng();
                edges.shuffle(&mut rng);
                let nremoved = nbonds.min(5).max(2);
                let n = rng.gen_range(1, nremoved);
                info!("will remove {} bonds ...", n);
                for i in 0..n {
                    mol.remove_bond(edges[i]);
                }
            } else {
                warn!("molecule has no bonds!");
            }

            let mol = kick(&mol)?;
            get_optimized_molecule(&mol, &self.config.runfile_opt)
        } else {
            Ok(mol.to_owned())
        }
    }
}

impl<'a> GeneticOperator for KickMutator<'a> {
    fn name() -> String {
        "Kick-start-Mutation".to_string()
    }
}

impl<'a> MutationOp<MoleculeGenome> for KickMutator<'a> {
    fn mutate<R>(&self, genome: MoleculeGenome, rng: &mut R) -> MoleculeGenome
    where
        R: Rng + Sized,
    {
        let mut mol = molecule_from_genome(&genome);
        let r: f64 = rng.gen_range(0.0, 1.0);
        let mol = if r < self.mutation_rate {
            info!("mutate molecule {:?}", mol.name);
            self.mutate_molecule(&mut mol, true)
                .map_err(|e| {
                    mol.to_file("/tmp/aa.mol2");
                    e
                })
                .expect("mutate molecule failed 1")
        } else {
            self.mutate_molecule(&mut mol, false)
                .map_err(|e| {
                    mol.to_file("/tmp/bb.mol2");
                    e
                })
                .expect("mutate molecule failed 2")
        };
        molecule_to_genome(&mol.sorted())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CutAndSpliceCrossBreeder<'a> {
    config: &'a Config,
}

impl<'a> GeneticOperator for CutAndSpliceCrossBreeder<'a> {
    fn name() -> String {
        "Cut-and-splice crossover".to_string()
    }
}

impl<'a> CrossoverOp<MoleculeGenome> for CutAndSpliceCrossBreeder<'a> {
    fn crossover<R>(
        &self,
        parents: Parents<MoleculeGenome>,
        rng: &mut R,
    ) -> Children<MoleculeGenome>
    where
        R: Rng + Sized,
    {
        info!("breed using crossover ...");
        let mol1 = molecule_from_genome(&parents[0]).sorted();
        let mol2 = molecule_from_genome(&parents[1]).sorted();
        let mut mol = spdkit::plane_cut_and_splice(&mol1, &mol2).expect("cut-and-splice failed");

        // avoid bad geometry which will cause opt failure
        mol.rebond();
        mol.clean().expect("clean");
        let mol = get_optimized_molecule(&mol, &self.config.runfile_opt)
            .map_err(|e| {
                error!("{:?}", e);
                e
            })
            .expect("crossover opt");
        let g = molecule_to_genome(&mol.sorted());

        vec![g]
    }
}
// operators:1 ends here

// genetic

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*genetic][genetic:1]]
// cluster structure search using genetic algorithm
pub fn genetic_search(config: &Config) -> Result<()> {
    let mrate = config.search.mutation_rate;
    info!("mutation rate: {}", mrate);

    let (initial_population, eref) = get_population_and_eref(&config);

    info!("reference energy = {}", eref);
    let mut feval = ClusterFitnessEvaluator::new(&config.runfile_sp, eref);
    feval.temperature = config.search.boltzmann_temperature;
    let algorithm = genetic_algorithm()
        .with_evaluation(feval.clone())
        // .with_selection(RouletteWheelSelector::new(0.7, 2))
        .with_selection(TournamentSelector::new(0.7, 2, 4, 1.0, true))
        .with_crossover(CutAndSpliceCrossBreeder {config: &config})
        .with_mutation(KickMutator {
            mutation_rate: mrate,
            config: &config,
        })
        .with_reinsertion(ElitistReinserter::new(
            feval,
            false,
            0.7,
        ))
        .with_initial_population(initial_population)
        .build();

    let ngen = config.search.max_generations;
    let mut magcalc_sim = simulate(algorithm)
        .until(GenerationLimit::new(ngen))
        .build();

    let mut i = 0;
    let mut old_fitness = 0;
    loop {
        match magcalc_sim.step() {
            Ok(SimResult::Intermediate(step)) => {
                let evaluated_population = step.result.evaluated_population;
                let average_fitness = evaluated_population.average_fitness();
                let best_solution = step.result.best_solution;
                let best_fitness = best_solution.solution.fitness;
                println!(
                    "Step: generation: {}, average_fitness: {}, \
                     best fitness: {},",
                    step.iteration, average_fitness, best_fitness
                );

                let mols: Vec<_> = evaluated_population
                    .individuals()
                    .iter()
                    .zip(evaluated_population.fitness_values())
                    .map(|(g, v)| molecule_from_genome(&g))
                    .collect();

                if *average_fitness > old_fitness {
                    old_fitness = *average_fitness;
                    i += 1;
                    let ofile = format!("./g{:03}.xyz", i);
                    io::write(ofile, &mols)?;
                }
            }
            Ok(SimResult::Final(step, processing_time, duration, stop_reason)) => {
                let best_solution = step.result.best_solution;
                println!("{}", stop_reason);
                println!(
                    "Final result after {}: generation: {}, \
                     best solution with fitness {} found in generation {}, processing_time: {}",
                    duration.fmt(),
                    step.iteration,
                    best_solution.solution.fitness,
                    best_solution.generation,
                    processing_time.fmt()
                );
                break;
            }
            Err(e) => {
                error!("{:?}", e);
                break;
            }
        }
    }

    Ok(())
}
// genetic:1 ends here
