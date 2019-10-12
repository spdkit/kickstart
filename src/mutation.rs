// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use gosh::gchemol;

use gchemol::prelude::*;
use gchemol::Molecule;
use spdkit::random::*;

use crate::common::*;
// imports:1 ends here

// random bond

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*random%20bond][random bond:1]]
pub fn random_bond_mutate(mol: &Molecule, degree: usize) -> Result<Molecule> {
    // import Educate trait here
    use educate::prelude::*;

    let mut new_mol = mol.clone();
    for _ in 0..degree {
        new_mol = educate::tmp_random_bond_mutate(&new_mol);
    }

    Ok(new_mol)
}
// random bond:1 ends here
