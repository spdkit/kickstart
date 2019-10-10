// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use crate::common::*;
use crate::config::Config;
use crate::core::*;
use crate::model::*;

use gchemol::prelude::*;
use gchemol::Molecule;
// imports:1 ends here

// initial seeds

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*initial%20seeds][initial seeds:1]]
/// Create a population with n individuals.
///
/// Initial molecule will be read from `molfile`
///
pub(crate) fn build_initial_genomes(config: &Config, n: Option<usize>) -> Vec<MolGenome> {
    info!("create initial population ..");
    let mol = Molecule::from_file(&config.molfile).expect("mol2");
    let n = n.unwrap_or(config.search.population_size);

    let mol_genomes: Vec<_> = (0..n)
        .into_par_iter()
        .map(|_| {
            crate::kick(&mol)
                .expect("kick mol")
                .get_optimized_molecule()
                .expect("opt mol")
                .encode()
        })
        .collect();

    mol_genomes
}
// initial seeds:1 ends here

// exploration: global search
// There are two approaches to create individuals in a gloabl sense:
// 1. cut-and-splice crossover
// 2. random-kick (initial-seeds)

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*exploration:%20global%20search][exploration: global search:1]]
pub(crate) fn global_add_new_genomes(n: usize) -> Vec<MolGenome> {
    // append randomly generated individuals
    let config = &crate::config::CONFIG;
    let random_gneomes = build_initial_genomes(&config, Some(n));
    random_gneomes
}

/// Create `n` molecules in random configurations using kickstart algorithm
///
/// # Parameters
///
/// * parent_mol: initial molecule containing multiple fragments (based on connectivity)
/// * n: the number of configurations to be generated
///
fn new_molecules_rand_kicked(parent_mol: &Molecule, n: usize) -> Vec<Molecule> {
    info!("Creating {} molecules using kickstart.", n);
    (0..n)
        .into_par_iter()
        .map(|_| crate::kick(&parent_mol).expect("kick parent_mol"))
        .collect()
}
// exploration: global search:1 ends here
