// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use crate::common::*;
use crate::core::*;
// imports:1 ends here

// public

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*public][public:1]]
pub(crate) fn new_random_genomes(n: usize) -> Vec<MolGenome> {
    use gchemol::prelude::*;
    use gchemol::Molecule;

    // FIXME: avoid global config below
    let config = &crate::config::CONFIG;
    let mol = Molecule::from_file(&config.molfile).expect("mol2");
    let mols = crate::kickstart::kick_bunch(&mol, n);

    // 1. optimize new molecule
    // 2. convert into genome
    crate::calculator::compute(mols)
        .expect("calc failure")
        .into_iter()
        .map(|mp| mp.encode())
        .collect()
}
// public:1 ends here
