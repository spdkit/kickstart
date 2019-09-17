// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use crate::common::*;
use crate::config::Config;

use gosh::gchemol::prelude::*;
use gosh::gchemol::{io, Atom, Molecule};

use spdkit::prelude::*;
// imports:1 ends here

// base

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*base][base:1]]
#[derive(Debug, Clone)]
struct MolIndividual;

/// The Genotype for molecule
#[derive(Clone, Debug)]
struct MolGenome {
    name: String,
    data: Vec<(usize, [f64; 3])>,
}

impl std::fmt::Display for MolGenome {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:}", self.name)
    }
}

impl spdkit::individual::Genome for MolGenome {}
// base:1 ends here

// bbm model

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*bbm%20model][bbm model:1]]
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
// bbm model:1 ends here

// evaluate

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*evaluate][evaluate:1]]
impl EvaluateObjectiveValue<MolGenome> for MolIndividual {
    fn evaluate(&self, genome: &MolGenome) -> f64 {
        let config = &crate::config::CONFIG;

        let mol = genome.decode();

        if let Ok(energy) = get_energy(&mol, &config.runfile_opt) {
            energy
        } else {
            warn!("no energy for {}", &mol.name);
            std::f64::MAX
        }
    }
}
// evaluate:1 ends here

// genome/molecule mapping
// genotype <=> phenotype conversion

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*genome/molecule%20mapping][genome/molecule mapping:1]]
trait ToGenome {
    fn encode(&self) -> MolGenome;
}

/// Create a random string of length `n` for naming a genome
fn random_name(n: usize) -> String {
    use rand::distributions::Alphanumeric;
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    rng.sample_iter(&Alphanumeric).take(n).collect()
}

impl ToGenome for Molecule {
    fn encode(&self) -> MolGenome {
        let mut g = vec![];
        for a in self.sorted().atoms() {
            let n = a.number();
            let p = a.position();
            g.push((n, p));
        }

        MolGenome {
            name: random_name(5),
            data: g,
        }
    }
}

impl MolGenome {
    fn decode(&self) -> Molecule {
        let mut mol = Molecule::new(&self.name);

        for (n, [x, y, z]) in &self.data {
            let a = Atom::build().element(*n).position(*x, *y, *z).finish();
            mol.add_atom(a);
        }

        mol.rebond();
        mol
    }
}
// genome/molecule mapping:1 ends here

// crossover

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*crossover][crossover:1]]
use spdkit::population::*;

#[derive(Debug, Clone)]
struct CutAndSpliceCrossOver;

impl VariationOperator<MolGenome> for CutAndSpliceCrossOver {
    fn breed_from<R: Rng + Sized>(
        &self,
        parents: &[Member<MolGenome>],
        rng: &mut R,
    ) -> Vec<MolGenome> {
        let config = &crate::config::CONFIG;

        let mol1 = parents[0].genome().decode();
        let mol2 = parents[1].genome().decode();
        info!("breed using crossover {} + {}.", mol1.name, mol2.name);

        let mut mol = crate::crossover::plane_cut_and_splice(&mol1, &mol2).expect("cut-and-splice failed");

        // avoid bad geometry which will cause opt failure
        mol.rebond();
        mol.clean().expect("clean");

        let mol = get_optimized_molecule(&mol, &config.runfile_opt)
            .with_context(|e| format!("failed with error: {:?}", e))
            .expect("crossover opt");

        vec![mol.encode()]
    }
}
// crossover:1 ends here

// mutation

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*mutation][mutation:1]]
impl Mutate for MolGenome {
    /// Mutate `n` bits randomly.
    fn mutate<R: Rng + Sized>(&mut self, n: usize, rng: &mut R) {
        let mutation_rate = crate::config::CONFIG.search.mutation_rate;

        let mut mol = self.decode();
        let r: f64 = rng.gen_range(0.0, 1.0);
        let mol = if r < mutation_rate {
            info!("mutate molecule {:?}", mol.name);
            mutate_molecule(&mut mol, true)
                .map_err(|e| {
                    mol.to_file("/tmp/aa.mol2");
                    e
                })
                .expect("mutate molecule failed 1")
        } else {
            mutate_molecule(&mut mol, false)
                .map_err(|e| {
                    mol.to_file("/tmp/bb.mol2");
                    e
                })
                .expect("mutate molecule failed 2")
        };

        *self = mol.encode();
    }
}

fn mutate_molecule(mol: &mut Molecule, mutate: bool) -> Result<Molecule> {
    let config = &crate::config::CONFIG;

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

        let mol = crate::kick(&mol)?;
        get_optimized_molecule(&mol, &config.runfile_opt)
    } else {
        Ok(mol.to_owned())
    }
}
// mutation:1 ends here

// initial seeds

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*initial%20seeds][initial seeds:1]]
/// create a population with n individuals.
/// initial molecule will be read from `molfile`
fn build_initial_genomes(config: &Config) -> Vec<MolGenome> {
    info!("create initial population ..");
    let mol = Molecule::from_file(&config.molfile).expect("mol2");
    let n = config.search.population_size;

    let mut mol_genomes = vec![];
    for i in 0..n {
        let mol = crate::kick(&mol).expect("kick mol");
        let mol = get_optimized_molecule(&mol, &config.runfile_opt).expect("opt mol");
        let e = get_energy(&mol, &config.runfile_sp).expect("sp energy");
        mol_genomes.push(mol.encode());
    }
    mol_genomes
}
// initial seeds:1 ends here

// genetic search

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*genetic%20search][genetic search:1]]
use spdkit::operators::selection::StochasticUniversalSampling as SusSelection;

// cluster structure search using genetic algorithm
pub fn genetic_search() -> Result<()> {
    let config = &crate::config::CONFIG;

    // create breeder gear
    let mrate = config.search.mutation_rate;
    info!("mutation rate: {}", mrate);
    let breeder = spdkit::GeneticBreeder::new()
        .with_selector(SusSelection::new(2))
        .with_crossover(CutAndSpliceCrossOver)
        .mutation_probability(mrate);

    // create valuer gear
    let temperature = config.search.boltzmann_temperature;
    let valuer = spdkit::Valuer::new()
        .with_fitness(spdkit::fitness::MinimizeEnergy::new(temperature))
        .with_creator(MolIndividual);

    // create evolution engine
    let mut engine = spdkit::Engine::new()
        .with_valuer(valuer)
        .with_breeder(breeder);

    // start the evolution loop
    let seeds = build_initial_genomes(&config);
    for g in engine.evolve(&seeds).take(config.search.max_generations) {
        let generation = g?;
        generation.summary();
        let energy = generation
            .population
            .best_member()
            .unwrap()
            .objective_value();
        if let Some(target_energy) = config.search.target_energy {
            if energy < target_energy {
                println!("target energy {} reached.", target_energy);
                break;
            }
        }

        // write generation results
        let mols: Vec<_> = generation
            .population
            .individuals()
            .iter()
            .map(|indv| indv.genome().decode())
            .collect();

        let ofile = format!("./g{:03}.xyz", generation.index);
        io::write(ofile, &mols)?;
    }

    Ok(())
}
// genetic search:1 ends here
