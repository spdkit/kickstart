// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use crate::common::*;
use crate::config::Config;
use crate::model::*;

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

// evaluate

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*evaluate][evaluate:1]]
impl EvaluateObjectiveValue<MolGenome> for MolIndividual {
    fn evaluate(&mut self, genome: &MolGenome) -> f64 {
        if let Ok(energy) = genome.decode().get_energy() {
            info!("evaluated indv {}, energy = {:-12.5}", genome.name, energy);
            energy
        } else {
            warn!("no energy for {}", genome.name);
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
        use educate::prelude::*;

        let mut mol = Molecule::new(&self.name);

        for (n, [x, y, z]) in &self.data {
            let a = Atom::build().element(*n).position(*x, *y, *z).finish();
            mol.add_atom(a);
        }

        mol.educated_rebond().unwrap();
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
        use educate::prelude::*;

        let config = &crate::config::CONFIG;

        let mol1 = parents[0].genome().decode();
        let mol2 = parents[1].genome().decode();
        debug!("breeding using crossover {} + {}.", mol1.name, mol2.name);

        let mut mol =
            crate::crossover::plane_cut_and_splice(&mol1, &mol2).expect("cut-and-splice failed");

        // avoid bad geometry which will cause opt failure
        // mol.rebond();
        mol.educate()
            .with_context(|_| format!("failed to educate"))
            .unwrap();

        let mol = mol.get_optimized_molecule().expect("crossover opt");

        let g = mol.encode();
        info!(
            "bred new indv {}, parents: {} + {}",
            g.name, mol1.name, mol2.name
        );

        vec![g]
    }
}
// crossover:1 ends here

// mutation

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*mutation][mutation:1]]
impl Mutate for MolGenome {
    /// Mutate `n` bits randomly.
    fn mutate<R: Rng + Sized>(&mut self, n: usize, rng: &mut R) {
        let mut mol = self.decode();
        mol = crate::mutation::random_bond_mutate(&mol)
            .expect("mutate molecule failed 1")
            .get_optimized_molecule()
            .expect("mutation opt");

        *self = mol.encode();
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
        let mol = mol.get_optimized_molecule().expect("opt mol");
        mol_genomes.push(mol.encode());
    }
    mol_genomes
}
// initial seeds:1 ends here

// breeder

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*breeder][breeder:1]]
#[derive(Clone)]
struct HyperMutation {
    mut_prob: f64,
}

impl Breed<MolGenome> for HyperMutation {
    /// Breed `m` new genomes from parent population.
    fn breed<R: Rng + Sized>(
        &mut self,
        m: usize,
        population: &Population<MolGenome>,
        rng: &mut R,
    ) -> Vec<MolGenome> {
        // loop until required number of genomes
        let mut required_genomes = Vec::with_capacity(m);
        while required_genomes.len() < m {
            // select 2 individuals as parents
            let selector = SusSelection::new(2);
            let parents = selector.select_from(population, rng);
            info!("selected {} parents for hyper mutation.", parents.len());

            for m in parents {
                let old_genome = m.genome();
                let old_score = m.objective_value();
                info!(
                    ">> Mutating parent {}, old_score = {}",
                    old_genome.name, old_score
                );

                // start mutation
                let mol = old_genome.decode();
                let mol = crate::mutation::random_bond_mutate(&mol)
                    .expect("mutate molecule failed")
                    .get_optimized_molecule()
                    .expect("mutation opt failed");

                required_genomes.push(mol.encode());
            }
        }

        required_genomes
    }
}
// breeder:1 ends here

// public

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*public][public:1]]
use spdkit::operators::selection::StochasticUniversalSampling as SusSelection;

// cluster structure search using genetic algorithm
pub fn genetic_search() -> Result<()> {
    let config = &crate::config::CONFIG;

    // create breeder gear
    let mrate = config.search.mutation_rate;
    info!("mutation rate: {}", mrate);
    // let breeder = spdkit::GeneticBreeder::new()
    //     .with_selector(SusSelection::new(2))
    //     .with_crossover(CutAndSpliceCrossOver)
    //     .mutation_probability(mrate);
    let breeder = HyperMutation { mut_prob: mrate };

    // create valuer gear
    let temperature = config.search.boltzmann_temperature;
    let valuer = spdkit::Valuer::new()
        .with_fitness(spdkit::fitness::MinimizeEnergy::new(temperature))
        .with_creator(MolIndividual);

    // create evolution engine
    let mut engine = spdkit::Engine::new()
        .with_valuer(valuer)
        .with_breeder(breeder);
    if let Some(n) = config.search.termination_nlast {
        println!("running mean termination: nlast = {}", n);
        engine.set_termination_nlast(n);
    };

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
// public:1 ends here
