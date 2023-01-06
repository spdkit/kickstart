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
    let config = &crate::config::CONFIG;
    for _ in 0..degree {
        if let Some(bbm_dir) = &config.md {
            let mut bbm = gosh::model::BlackBoxModel::from_dir(bbm_dir)?;
            info!("Exploit using minima hopping algorithm with bbm: {bbm_dir:?}");
            // FIXME: removed
            // mutate_by_md(&mut new_mol, &mut bbm)?;
        }
        new_mol = educate::adhoc::random_bond_mutate(&new_mol);
    }

    Ok(new_mol)
}
// 4ebcd3f9 ends here
