// [[file:../kickstart.note::*imports][imports:1]]
use gosh::gchemol;

use gchemol::prelude::*;
use gchemol::Molecule;
use spdkit::random::*;

use crate::common::*;
// imports:1 ends here

// [[file:../kickstart.note::4ebcd3f9][4ebcd3f9]]
pub fn random_bond_mutate(mol: &Molecule, degree: usize) -> Result<Molecule> {
    // import Educate trait here
    use educate::prelude::*;

    let mut new_mol = mol.clone();
    for _ in 0..degree {
        new_mol = educate::adhoc::random_bond_mutate(&new_mol);
    }

    Ok(new_mol)
}
// 4ebcd3f9 ends here
