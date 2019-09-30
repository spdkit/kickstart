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
// import Educate trait here
use educate::prelude::*;

pub fn random_bond_mutate(mol: &Molecule, degree: usize) -> Result<Molecule> {
    use crate::model::*;

    let mut new_mol = mol.clone();
    for _ in 0..degree {
        new_mol = educate::tmp_random_bond_mutate(&new_mol);
    }

    Ok(new_mol)
}
// random bond:1 ends here

pub(crate) fn mutate_molecule(mol: &Molecule) -> Result<Molecule> {
    let mut mol = mol.clone();
    let nbonds = mol.nbonds();
    if nbonds > 0 {
        // remove bonds randomly to break molecule into parts
        let mut edges: Vec<_> = mol.bonds().map(|b| b.index()).collect();
        let mut rng = thread_rng();
        edges.shuffle(&mut rng);
        let nremoved = nbonds.min(5).max(2);
        let n = rng.gen_range(1, nremoved);
        info!("will remove {} bonds ...", n);
        for i in 0..n {
            mol.remove_bond(edges[i]);
        }
    } else {
        warn!("molecule has no bonds!");
    }

    crate::kick(&mol)
}
