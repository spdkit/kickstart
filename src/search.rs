// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use std::iter::FromIterator;

use genevo::operator::prelude::*;
use genevo::prelude::*;
use genevo::types::fmt::Display;

use gosh::gchemol::prelude::*;
use gosh::gchemol::{Atom, Molecule};

use crate::common::*;
use crate::config::Config;
// imports:1 ends here

// base

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*base][base:1]]
// energy => fitness
fn calc_fitness(energy: f64) -> u32 {
    let temperature = 10000.;
    let eref = -13.0;
    let value = (energy - eref) * 96.;
    let fitness = (-1.0 * value / (temperature * 0.0083145)).exp();

    (fitness * 1000.) as u32
}

#[test]
fn test_calc_fit() {
    let x = -205.30249;
    let x = calc_fitness(x);
    println!("{:#?}", x);
}

type MoleculeGenome = Vec<(usize, f64, f64, f64)>;

fn molecule_from_genome(genome: &MoleculeGenome) -> Molecule {
    let mut mol = Molecule::new("genome");
    for (n, x, y, z) in genome {
        let a = Atom::build().element(*n).position(*x, *y, *z).finish();
        mol.add_atom(a);
    }

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

// fitness

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*fitness][fitness:1]]
use gosh::models::*;

#[derive(Clone, Debug, PartialEq)]
struct ClusterFitnessEvaluator<'a> {
    runfile: &'a str,
}

impl<'a> ClusterFitnessEvaluator<'a> {
    fn new(runfile: &'a str) -> Self {
        ClusterFitnessEvaluator { runfile: &runfile }
    }
}

// calculate geometry and energy using bbm model
fn get_energy(mol: &Molecule, runfile: &str) -> Result<f64> {
    info!("calculate single point energy using bbm model ...");

    let mut bbm = BlackBox::from_dir(runfile);
    let mr = bbm.compute(mol)?;

    if let Some(energy) = mr.energy {
        println!("final energy = {:-12.5}", energy);
        return Ok(energy);
    } else {
        bail!("no energy record found in the output!");
    }
}

fn get_optimized_molecule(mol: &Molecule, runfile: &str) -> Result<Molecule> {
    info!("optimize geometry using bbm model ...");
    let mut bbm = BlackBox::from_dir(runfile);
    dbg!(runfile);
    let mr = bbm.compute(mol)?;

    if let Some(energy) = mr.energy {
        println!("opt energy = {:-12.5}", energy);

        info!("output final molecule ...");
        if let Some(mut mol) = mr.molecule {
            mol.rebond();
            return Ok(mol);
        } else {
            bail!("no molecule record found in bbm results");
        }
    } else {
        bail!("no energy record found in the output!");
    }
}

impl<'a> FitnessFunction<MoleculeGenome, u32> for ClusterFitnessEvaluator<'a> {
    fn fitness_of(&self, individual: &MoleculeGenome) -> u32 {
        info!("evaluate fitness ...");
        let mol = molecule_from_genome(&individual);
        // mapping energy to fitness
        let fitness = {
            if let Ok(energy) = get_energy(&mol, &self.runfile) {
                calc_fitness(energy)
            } else {
                warn!("no energy for {}", &mol.name);
                self.highest_possible_fitness()
            }
        };

        fitness
    }

    fn average(&self, fitness_values: &[u32]) -> u32 {
        (fitness_values.iter().sum::<u32>() / fitness_values.len() as u32)
    }

    fn highest_possible_fitness(&self) -> u32 {
        1000
    }

    fn lowest_possible_fitness(&self) -> u32 {
        0
    }
}
// fitness:1 ends here

// operators

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*operators][operators:1]]
use genevo::genetic::{Children, Genotype, Parents};
use genevo::operator::prelude::*;
use genevo::operator::{CrossoverOp, GeneticOperator, MutationOp};
use genevo::prelude::*;
use genevo::random::Rng;
use rand::prelude::*;
use spdkit;

use crate::kick;

// create a population with n individuals
// initial molecule will be read from `molfile`
fn create_population(config: &Config) -> Population<MoleculeGenome> {
    info!("create population using kick ..");
    let mol = Molecule::from_file(&config.molfile).expect("mol2");
    let mut individuals = vec![];
    for i in 0..10 {
        let mol = kick(&mol).expect("kick mol");
        dbg!(i);
        let mol = get_optimized_molecule(&mol, &config.runfile_opt).expect("opt mol");
        let g = molecule_to_genome(&mol);
        individuals.push(g);
    }

    Population::with_individuals(individuals)
}

#[derive(Clone, Debug, PartialEq)]
pub struct KickMutator<'a> {
    mutation_rate: f64,
    config: &'a Config,
}

impl<'a> KickMutator<'a> {
    pub fn mutate_molecule(&self, mol: &mut Molecule, mutate: bool) -> Result<Molecule> {
        let nbonds = mol.nbonds();
        info!("mutate molecule {:?}", mol.name);

        if nbonds > 0 {
            // choose n bonds randomly
            let mut rng = thread_rng();
            let mut n = rng.gen_range(0, nbonds);

            let edges: Vec<_> = mol.bonds().map(|b| b.index()).collect();

            mol.remove_bond(edges[n]);
        } else {
            info!("molecule has no bonds.");
        }

        if mutate {
            let mol = kick(&mol)?;
            get_optimized_molecule(&mol, &self.config.runfile_opt)
        } else {
            get_optimized_molecule(&mol, &self.config.runfile_opt)
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
            self.mutate_molecule(&mut mol, true)
                .expect("mutate molecule failed")
        } else {
            self.mutate_molecule(&mut mol, false)
                .expect("mutate molecule failed")
        };
        molecule_to_genome(&mol)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CutAndSpliceCrossBreeder {}

impl GeneticOperator for CutAndSpliceCrossBreeder {
    fn name() -> String {
        "Cut-and-splice crossover".to_string()
    }
}

impl CrossoverOp<MoleculeGenome> for CutAndSpliceCrossBreeder {
    fn crossover<R>(
        &self,
        parents: Parents<MoleculeGenome>,
        rng: &mut R,
    ) -> Children<MoleculeGenome>
    where
        R: Rng + Sized,
    {
        let mut mol1 = molecule_from_genome(&parents[0]).sorted();
        let mut mol2 = molecule_from_genome(&parents[1]).sorted();
        let mol = spdkit::plane_cut_and_splice(&mol1, &mol2).expect("cut-and-splice failed");

        let g = molecule_to_genome(&mol.sorted());

        vec![g]
    }
}
// operators:1 ends here

// main loop

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*main%20loop][main loop:1]]
// cluster structure search using genetic algorithm
pub fn genetic_search(config: &Config) -> Result<()> {
    let initial_population = create_population(&config);

    // runfile for final energy calculation
    let runfile = &config.runfile_sp;

    let algorithm = genetic_algorithm()
        .with_evaluation(ClusterFitnessEvaluator::new(&runfile))
        .with_selection(RouletteWheelSelector::new(0.7, 2))
        .with_crossover(CutAndSpliceCrossBreeder {})
        .with_mutation(KickMutator {
            mutation_rate: 0.1,
            config: &config,
        })
        .with_reinsertion(ElitistReinserter::new(
            ClusterFitnessEvaluator::new(&runfile),
            false,
            0.7,
        ))
        .with_initial_population(initial_population)
        .build();

    let ngen = &config.search.max_generations;
    let mut magcalc_sim = simulate(algorithm)
        .until(GenerationLimit::new(*ngen))
        .build();

    let mut igen = 1;
    loop {
        println!("generation: {:}", igen);
        let result = magcalc_sim.step();
        match result {
            Ok(SimResult::Intermediate(step)) => {
                let evaluated_population = step.result.evaluated_population;
                let best_solution = step.result.best_solution;
                println!(
                    "Step: generation: {}, average_fitness: {}, \
                     best fitness: {}, duration: {}, processing_time: {}",
                    step.iteration,
                    evaluated_population.average_fitness(),
                    best_solution.solution.fitness,
                    step.duration.fmt(),
                    step.processing_time.fmt()
                );
                println!("{:#?}", best_solution);
                let g = best_solution.solution.genome;
                let mut mol = molecule_from_genome(&g);
                let ofile = format!("./best-g{:03}.mol2", step.iteration);
                mol.rebond();
                mol.to_file(&ofile);
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
            Err(error) => {
                println!("{}", error.display());
                break;
            }
        }
        igen += 1;
    }

    Ok(())
}
// main loop:1 ends here
