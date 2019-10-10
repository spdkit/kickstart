// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use crate::common::*;
use crate::config::Config;
use crate::core::*;
use crate::exploitation::*;
use crate::exploration::*;
use crate::model::*;

use gosh::gchemol::prelude::*;
use gosh::gchemol::{io, Atom, Molecule};

use spdkit::operators::selection::StochasticUniversalSampling as SusSelection;
use spdkit::prelude::*;
use spdkit::*;
// imports:1 ends here

// evaluate

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*evaluate][evaluate:1]]
impl EvaluateObjectiveValue<MolGenome> for MolIndividual {
    fn evaluate(&self, genome: &MolGenome) -> f64 {
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

// evolve
// core part starts here

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
    let p_add_global_kick = 0.2;
    let p_add_global_crossover = 0.3;

    // let n_local_search = (n as f64 * 0.8) as usize;
    let n_local_search = n;
    info!("vds: select {} genomes for exploitation", n_local_search);
    let selector = SusSelection::new(2);
    let mut next_parent = || {
        let mut rng = rand::thread_rng();
        // add new genome locally
        if rng.gen::<f64>() < p_add_global_kick {
            let new_genome = global_add_new_genomes(1).pop().unwrap();
            info!("candidate genome {}: randomly generated", new_genome.name);
            new_genome
        }
        // add new genome globally
        else {
            let members = selector.select_from(cur_population, &mut rng);
            let new_genome =
            // global add using cut-and-splice crossover
            if rng.gen::<f64>() < p_add_global_crossover {
                let new_genome = CutAndSpliceCrossOver
                    .breed_from(&members, &mut rng)
                    .pop()
                    .unwrap();
                info!("candidate genome {}: cut-and-splice crossover", new_genome.name);
                new_genome
            }
            // global add using random kick
            else {
                let new_genome = members[0].genome().to_owned();
                info!("candidate genome {}: weighted selection", new_genome.name);
                new_genome
            };
            let members = selector.select_from(cur_population, &mut rng);
            new_genome
        }
    };

    // vds exploitation in parallel
    let required_genomes: Vec<_> = (0..n_local_search)
        .into_par_iter()
        .map(|_| {
            let parent = next_parent();
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
        let new_genomes = global_add_new_genomes(n_add_random);
        let new_indvs = valuer.create_individuals(new_genomes);
        all_indvs.extend_from_slice(&new_indvs);
    }

    valuer.build_population(all_indvs[..m].to_vec())
}
// evolve:1 ends here

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
// public:1 ends here
