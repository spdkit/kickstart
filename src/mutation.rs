// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use gosh::gchemol;

use gchemol::prelude::*;
use gchemol::Molecule;
use spdkit::random::*;

use crate::common::*;
// imports:1 ends here

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

// random bond

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*random%20bond][random bond:1]]
// import Educate trait here
use educate::prelude::*;

pub(crate) fn random_bond_mutate(mol: &Molecule) -> Result<Molecule> {
    use gchemol::Bond;

    // randomly bond a pair of atoms
    let nodes: Vec<_> = mol.atoms().map(|a| a.index()).collect();

    // choose two different nodes
    let mut rng = thread_rng();
    let mut random_choose_node = || nodes.choose(&mut rng).expect("no node");
    let node1 = random_choose_node();
    let connected: Vec<_> = mol.neighbors(*node1);

    // avoid infinite loop
    let mut node2 = random_choose_node();
    assert!(nodes.len() >= 2);
    // exclude current node and its neighboring nodes
    while node1 == node2 || connected.contains(node2) {
        node2 = random_choose_node();
    }

    let mut mol = mol.clone();
    mol.add_bond(*node1, *node2, Bond::single());
    // dbg!(node1, node2);
    mol.educated_clean();

    Ok(mol)
}
// random bond:1 ends here

// test

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*test][test:1]]
#[test]
fn test_rand_bond_mutate() -> Result<()> {
    let mol = Molecule::from_file("/tmp/test1.mol2")?;

    let mols: Vec<_> = (0..10).map(|_| random_bond_mutate(&mol).unwrap()).collect();

    gchemol::io::write("/tmp/a.mol2", &mols)?;
    Ok(())
}
// test:1 ends here
