// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use gosh::gchemol;

use gchemol::prelude::*;
use gchemol::Molecule;
use spdkit::random::*;

use crate::common::*;
// imports:1 ends here

// random bond / new

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*random%20bond%20/%20new][random bond / new:1]]
// import Educate trait here
use educate::prelude::*;

pub(crate) fn random_bond_mutate(mol: &Molecule) -> Result<Molecule> {
    let mut mol = educate::tmp_random_bond_mutate(mol);

    mol.educated_clean();

    Ok(mol)
}
// random bond / new:1 ends here

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

// test

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*test][test:1]]
#[test]
#[ignore]
fn test_rand_bond_mutate() -> Result<()> {
    let mol = Molecule::from_file("/tmp/test1.mol2")?;

    let mols: Vec<_> = (0..10).map(|_| random_bond_mutate(&mol).unwrap()).collect();

    gchemol::io::write("/tmp/a.mol2", &mols)?;
    Ok(())
}
// test:1 ends here
