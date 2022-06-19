// [[file:../kickstart.note::*imports][imports:1]]
use crate::common::*;
use crate::core::*;
// imports:1 ends here

// [[file:../kickstart.note::d22e2006][d22e2006]]
pub fn new_random_genomes(n: usize) -> Vec<MolGenome> {
    info!("create {n} random genomes ...");
    use gchemol::prelude::*;
    use gchemol::Molecule;

    // FIXME: avoid global config below
    let config = &crate::config::CONFIG;
    let mol = Molecule::from_file(&config.molfile).expect("mol2");
    let mols = crate::kickstart::kick_bunch(&mol, n);
    info!("Created {} random molecules.", mols.len());

    crate::model::evaluate_molecules(mols)
        .expect("comput failure")
        .into_iter()
        .map(|e| e.genome)
        .collect()
}
// d22e2006 ends here
