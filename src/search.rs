// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use crate::common::*;
use crate::config::Config;
use crate::core::*;
use crate::exploitation::*;
use crate::model::*;

use gosh::gchemol::prelude::*;
use gosh::gchemol::{io, Atom, Molecule};

use spdkit::operators::selection::StochasticUniversalSampling as SusSelection;
use spdkit::prelude::*;
use spdkit::*;
// imports:1 ends here

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

        let mol1 = parents[0].genome().decode();
        let mol2 = parents[1].genome().decode();
        debug!("breeding using crossover {} + {}.", mol1.name, mol2.name);

        let mut mol =
            crate::crossover::plane_cut_and_splice(&mol1, &mol2).expect("cut-and-splice failed");

        // avoid bad geometry which will cause opt failure
        mol.educate()
            .with_context(|_| format!("failed to educate"))
            .unwrap();

        let genomes: Vec<_> = crate::calculator::compute(vec![mol])
            .expect("calc failure")
            .into_iter()
            .map(|mp| mp.encode())
            .collect();

        info!(
            "bred {}, parents: {} + {}",
            &genomes[0], mol1.name, mol2.name
        );

        genomes
    }
}
// crossover:1 ends here

// next parent

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*next%20parent][next parent:1]]
fn next_parent(cur_population: &Population<MolGenome>) -> MolGenome {
    let p_add_global_kick = 0.2;
    let p_add_global_crossover = 0.3;

    let mut rng = rand::thread_rng();
    let selector = SusSelection::new(2);

    // add new genome globally
    if rng.gen::<f64>() < p_add_global_kick {
        let new_genome = crate::exploration::new_random_genomes(1).pop().unwrap();
        info!("candidate genome {}: randomly generated", new_genome);
        new_genome
    }
    // add new genome globally based on crossover
    else {
        let members = selector.select_from(cur_population, &mut rng);
        let new_genome =
            // global add using cut-and-splice crossover
            if rng.gen::<f64>() < p_add_global_crossover {
                let new_genome = CutAndSpliceCrossOver
                    .breed_from(&members, &mut rng)
                    .pop()
                    .unwrap();
                info!("candidate genome {}: cut-and-splice crossover", new_genome);
                new_genome
            }
            // add new genome locally using random kick
            else {
                let new_genome = members[0].genome().to_owned();
                info!("candidate genome {}: weighted selection", new_genome);
                new_genome
            };
        info!("candidate genome {}: weighted selection", new_genome);
        new_genome
    }
}
// next parent:1 ends here

// evolve
// core part starts here

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*evolve][evolve:1]]
struct MyAlgorithm {}

impl Default for MyAlgorithm {
    fn default() -> Self {
        Self {}
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
    let similarity_energy_threshold = 0.01; // in eV, adhoc
    let m = cur_population.size_limit();
    let n = cur_population.size();
    let p_add_global_kick = 0.2;
    let p_add_global_crossover = 0.3;

    let n_local_search = n;
    info!("vds: select {} genomes for exploitation", n_local_search);
    // vds exploitation in parallel
    let required_genomes: Vec<_> = (0..n_local_search)
        .into_par_iter()
        .map(|_| {
            let parent = next_parent(cur_population);
            variable_depth_search(&parent)
        })
        .collect();

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
        let new_genomes = crate::exploration::new_random_genomes(n_add_random);
        let new_indvs = valuer.create_individuals(new_genomes);
        all_indvs.extend_from_slice(&new_indvs);
    }

    valuer.build_population(all_indvs[..m].to_vec())
}
// evolve:1 ends here

// survive

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*survive][survive:1]]
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
    let algo = MyAlgorithm::default();

    // create evolution engine
    let mut engine = spdkit::Engine::create().valuer(valuer).algorithm(algo);
    if let Some(n) = config.search.termination_nlast {
        println!("running mean termination: nlast = {}", n);
        engine.set_termination_nlast(n);
    };

    // start the evolution loop
    let nbunch = config.search.population_size;
    let seeds = crate::exploration::new_random_genomes(nbunch);
    for g in engine.evolve(&seeds).take(config.search.max_generations) {
        let generation = g?;
        generation.summary();
        let energy = generation
            .population
            .best_member()
            .unwrap()
            .objective_value();

        // write generation results
        let mols: Vec<_> = generation
            .population
            .individuals()
            .iter()
            .map(|indv| indv.genome().decode())
            .collect();

        let ofile = format!("./g{:03}.xyz", generation.index);
        io::write(ofile, &mols)?;

        if let Some(target_energy) = config.search.target_energy {
            if energy < target_energy {
                println!("target energy {} reached.", target_energy);
                break;
            }
        }
    }

    Ok(())
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
// public:1 ends here
