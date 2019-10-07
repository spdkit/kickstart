// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use crate::common::*;
use crate::config::Config;
use crate::model::*;

use gosh::gchemol::prelude::*;
use gosh::gchemol::{io, Atom, Molecule};

use spdkit::prelude::*;
use spdkit::*;
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
        let mol = self.decode();
        let mol = crate::mutation::random_bond_mutate(&mol, n)
            .expect("mutate molecule failed")
            .get_optimized_molecule()
            .expect("mutation opt failed");
        *self = mol.encode();
    }
}
// mutation:1 ends here

// initial seeds

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*initial%20seeds][initial seeds:1]]
/// create a population with n individuals.
/// initial molecule will be read from `molfile`
fn build_initial_genomes(config: &Config, n: Option<usize>) -> Vec<MolGenome> {
    info!("create initial population ..");
    let mol = Molecule::from_file(&config.molfile).expect("mol2");
    let n = n.unwrap_or(config.search.population_size);

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
use spdkit::operators::selection::StochasticUniversalSampling as SusSelection;

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
                let mol = crate::mutation::random_bond_mutate(&mol, 1)
                    .expect("mutate molecule failed")
                    .get_optimized_molecule()
                    .expect("mutation opt failed");

                required_genomes.push(mol.encode());
            }
        }

        // append randomly generated individuals
        let config = &crate::config::CONFIG;
        let random_gneomes = build_initial_genomes(&config, Some(5));
        required_genomes.extend_from_slice(&random_gneomes);

        required_genomes
    }
}
// breeder:1 ends here

// evolve

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*evolve][evolve:1]]
struct MyAlgorithm {
    //
}

impl Default for MyAlgorithm {
    fn default() -> Self {
        Self {
            // --
        }
    }
}

impl<F, C> Evolve<MolGenome, F, C> for MyAlgorithm
where
    F: EvaluateFitness<MolGenome>,
    C: EvaluateObjectiveValue<MolGenome>,
{
    fn next_generation(
        &mut self,
        cur_population: &Population<MolGenome>,
        valuer: &mut Valuer<MolGenome, F, C>,
    ) -> Population<MolGenome> {
        evolve_core(cur_population, valuer)
    }
}

fn evolve_core<F, C>(
    cur_population: &Population<MolGenome>,
    valuer: &mut Valuer<MolGenome, F, C>,
) -> Population<MolGenome>
where
    F: EvaluateFitness<MolGenome>,
    C: EvaluateObjectiveValue<MolGenome>,
{
    // constants first
    let similarity_energy_threshold = 0.01; // in eV
    let m = cur_population.size_limit();
    let n = cur_population.size();
    let p_add_global = 0.2;
    let n_local_search = (n as f64 * 0.4) as usize;

    info!("vds: select {} genomes for exploitation", n_local_search);
    let mut required_genomes = vec![];
    let mut rng = rand::thread_rng();
    let selector = SusSelection::new(1);
    while required_genomes.len() < n_local_search {
        // add one randomly generated genome
        let parent = if rng.gen::<f64>() < p_add_global {
            let new_genome = global_add_new_genomes(1).pop().unwrap();
            info!("candidate genome {}: randomly generated", new_genome.name);
            new_genome
        } else {
            let members = selector.select_from(cur_population, &mut rng);
            let new_genome = members[0].genome().to_owned();
            info!("candidate genome {}: weighted selection", new_genome.name);
            new_genome
        };
        let new_genome = variable_depth_search(&parent);
        required_genomes.push(new_genome);
    }

    // create a new population from old population and new generated genomes
    let old_genomes: Vec<_> = cur_population
        .individuals()
        .iter()
        .map(|indv| indv.genome().to_owned())
        .collect();
    let all_genomes = [required_genomes, old_genomes].concat();
    let mut all_indvs = valuer.create_individuals(all_genomes);
    // remove similar individuals
    all_indvs.remove_duplicates_by_energy(similarity_energy_threshold);

    // Add new random genomes only when it really needs. This will reduce
    // redudant calculations.
    if all_indvs.len() < m {
        let n_add_random = m - all_indvs.len();
        // add random genomes as candicates in a global way
        println!("Add {} new indvs with random kick", n_add_random);
        let new_genomes = global_add_new_genomes(n_add_random);
        let new_indvs = valuer.create_individuals(new_genomes);
        all_indvs.extend_from_slice(&new_indvs);
    }

    valuer.build_population(all_indvs[..m].to_vec())
}
// evolve:1 ends here

// global search
// There are two approaches to create individuals in a gloabl sense:
// 1. cut-and-splice crossover
// 2. random-kick (initial-seeds)

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*global%20search][global search:1]]
fn global_add_new_genomes(n: usize) -> Vec<MolGenome> {
    // append randomly generated individuals
    let config = &crate::config::CONFIG;
    let random_gneomes = build_initial_genomes(&config, Some(n));
    random_gneomes
}
// global search:1 ends here

// local search: variable depth search

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*local%20search:%20variable%20depth%20search][local search: variable depth search:1]]
/// Return the best MolGenome during local search neighbors from `old_genome`.
fn variable_depth_search(old_genome: &MolGenome) -> MolGenome {
    use crate::model::*;
    use spdkit::common::float_ordering_minimize;

    // search options
    let search_width = 3;
    let search_depth = 3;
    let max_iterations = 6;

    // core iteration for vds
    let iteration = || {
        // start mutation
        let mut mutated: Vec<_> = (0..search_width)
            .map(|_| old_genome.mutated_with_energy())
            .collect();
        // take current best
        mutated.sort_by(|a, b| float_ordering_minimize(&a.1, &b.1));
        let best = &mutated[0].0;
        mutated[0].clone()
    };

    let ini_energy = old_genome.decode().energy();
    let mut old_energy = ini_energy;
    info!("initial genome {} energy = {}", old_genome.name, old_energy);
    let mut search_results = vec![(old_genome.clone(), old_energy)];
    let mut istop = 0;
    for icycle in 0.. {
        info!("{:=^72}", format!(" vds iteration {} ", icycle));
        let (new_genome, new_energy) = iteration();
        info!(
            "current best in {} candidates: {} energy = {}",
            search_width, new_genome.name, new_energy
        );

        if new_energy > old_energy {
            istop += 1;
        }
        // reset stop flag when energy is decreasing
        else {
            istop = 0;
        }

        // stop iterations when situation get worse
        if istop >= search_depth {
            break;
        }
        old_energy = new_energy;
        search_results.push((new_genome, new_energy));
        if icycle >= max_iterations {
            info!("Max allowed iteration reached during vds.");
            break;
        }
    }

    // select the best one during search.
    search_results.sort_by(|a, b| float_ordering_minimize(&a.1, &b.1));

    let best = search_results.remove(0);
    println!(
        "vds result: {} => {}, {:-10.4} => {:-10.4} ",
        old_genome.name, best.0.name, ini_energy, best.1,
    );
    best.0
}

impl MolGenome {
    /// Return a mutated MolGenome using rand-bond-mutation algorithm.
    fn mutated_with_energy(&self) -> (Self, f64) {
        let mol = self.decode();
        let mol = crate::mutation::random_bond_mutate(&mol, 1)
            .expect("mutate molecule failed")
            .get_optimized_molecule()
            .expect("mutation opt failed");

        let energy = mol.energy();
        let genome = mol.encode();
        debug!("mutant energy of {} = {}", genome.name, energy);
        (genome, energy)
    }
}

#[test]
#[ignore]
fn test_vds() -> Result<()> {
    use gchemol::prelude::*;
    use gosh::gchemol;

    let mols = gchemol::io::read("/tmp/g000.xyz")?;
    let ini_genome = mols[0].encode();
    let _ = variable_depth_search(&ini_genome);
    Ok(())
}
// local search: variable depth search:1 ends here

// survive

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*survive][survive:1]]
use spdkit::prelude::*;
use spdkit::{Genome, Individual, Population};

#[derive(Clone)]
pub struct Survivor;

impl Survive<MolGenome> for Survivor {
    fn survive<R: Rng + Sized>(
        &mut self,
        population: Population<MolGenome>,
        _rng: &mut R,
    ) -> Vec<Individual<MolGenome>> {
        get_survived_individuals(population)
    }
}

/// remove dumplicates and weak individuals to fit population size limit.
fn get_survived_individuals(population: Population<MolGenome>) -> Vec<Individual<MolGenome>> {
    // FIXME: adhoc hacking for removing duplicates, based on energy
    // criterion only
    let threshold = 0.01;

    let n_old = population.size();
    let mut members: Vec<_> = population.members().collect();
    members.sort_by_fitness();

    // FIXME: adhoc
    let mut to_keep: Vec<_> = members.into_iter().enumerate().collect();
    to_keep.dedup_by(|a, b| {
        let (ma, mb) = (&a.1, &b.1);
        let (va, vb) = (ma.objective_value(), mb.objective_value());
        (va - vb).abs() < threshold
    });
    let n_remove = n_old - to_keep.len();
    info!("removed {} duplicates", n_remove);

    let mut indvs = vec![];
    for p in to_keep.into_iter().take(population.size_limit()) {
        let m = p.1;
        indvs.push(m.individual.to_owned());
    }

    indvs
}

trait RemoveDuplicates {
    fn remove_duplicates_by_energy(&mut self, threshold: f64) -> usize;
}

// FIXME: adhoc hacking for removing duplicates, based on energy
// criterion only
impl RemoveDuplicates for Vec<Individual<MolGenome>> {
    /// # Parameters
    /// * threshold: energy diff threshold for similarity detection
    ///
    /// # Returns
    /// * return the number of removed indvs
    fn remove_duplicates_by_energy(&mut self, threshold: f64) -> usize {
        use spdkit::common::float_ordering_minimize;

        let n_old = self.len();
        self.sort_by(|a, b| float_ordering_minimize(&a.objective_value(), &b.objective_value()));

        // FIXME: adhoc
        self.dedup_by(|a, b| (a.objective_value() - b.objective_value()).abs() < threshold);
        let n_removed = n_old - self.len();
        info!("Removed {} duplicates", n_removed);
        n_removed
    }
}
// survive:1 ends here

// public

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*public][public:1]]
// cluster structure search using genetic algorithm
pub fn genetic_search() -> Result<()> {
    let config = &crate::config::CONFIG;

    // create valuer gear
    let temperature = config.search.boltzmann_temperature;
    let valuer = spdkit::Valuer::new()
        .with_fitness(spdkit::fitness::MinimizeEnergy::new(temperature))
        .with_creator(MolIndividual);

    // setup evolution algorithm
    // create breeder gear
    // let mrate = config.search.mutation_rate;
    // info!("mutation rate: {}", mrate);
    // let breeder = HyperMutation { mut_prob: mrate };
    // // create a survivor gear
    // let survivor = Survivor;
    // let algo = spdkit::EvolutionAlgorithm::new(breeder, survivor);
    let algo = MyAlgorithm::default();

    // create evolution engine
    let mut engine = spdkit::Engine::create().valuer(valuer).algorithm(algo);

    if let Some(n) = config.search.termination_nlast {
        println!("running mean termination: nlast = {}", n);
        engine.set_termination_nlast(n);
    };

    // start the evolution loop
    let seeds = build_initial_genomes(&config, None);
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
