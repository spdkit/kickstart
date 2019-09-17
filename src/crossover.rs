// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use std::collections::HashSet;

use gosh::gchemol;

use gchemol::prelude::*;
use gchemol::Molecule;
use gchemol::{geometry::rand_rotate};

use crate::common::*;
// imports:1 ends here

// core

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*core][core:1]]
// return indices of atoms lying above the cutting plane
fn indices_above_plane(positions: &Vec<[f64; 3]>) -> HashSet<usize> {
    positions
        .iter()
        .enumerate()
        .filter(|(i, p)| p[2].is_sign_positive())
        .map(|(i, p)| i)
        .collect()
}

// cut the molecule into two parts using a random plane
fn cut_molecule_by_rand_plane(mol: &Molecule) -> (HashSet<usize>, HashSet<usize>, Vec<[f64; 3]>) {
    let natoms = mol.natoms();

    let mut mol = mol.clone();
    mol.recenter();

    let positions = mol.positions();
    let rotated = rand_rotate(&positions);

    let ind_all: HashSet<_> = (0..natoms).collect();
    let ind_above = indices_above_plane(&rotated);
    let ind_below = ind_all.difference(&ind_above).map(|x| *x).collect();

    (ind_above, ind_below, rotated)
}

pub(crate) fn plane_cut_and_splice(mol1: &Molecule, mol2: &Molecule) -> Result<Molecule> {
    let natoms = mol1.natoms();
    // sanity check
    if mol2.natoms() == natoms {
        for (a1, a2) in mol1.atoms().zip(mol2.atoms()) {
            if a1.symbol() != a2.symbol() {
                bail!("atom numbering is inconsistent!");
            }
        }
    } else {
        bail!("molecules in difference size.");
    }

    // record element symbols
    let symbols = mol1.symbols();
    let reduced_symbols = mol1.reduced_symbols();
    let maxloop = 50000;
    let mut iloop = 0;

    let mut omol = mol1.clone();
    loop {
        let (above1, below1, rotated1) = cut_molecule_by_rand_plane(&mol1);
        let (above2, below2, rotated2) = cut_molecule_by_rand_plane(&mol2);
        debug!("above1 {:?}", above1);
        debug!("below2 {:?}", below2);
        // check if number of atoms is correct
        if above1.len() + below2.len() == natoms {
            // check if element types is correct

            let s1: Vec<_> = above1.iter().map(|&i| symbols[i]).collect();

            let mut s2: Vec<_> = below2.iter().map(|&i| symbols[i]).collect();

            s2.extend(s1.iter());

            assert_eq!(natoms, s2.len());
            omol.set_symbols(&s2);
            let mut got = true;
            for (k, v) in omol.reduced_symbols() {
                if reduced_symbols[&k] != v {
                    got = false;
                }
            }
            if got {
                // update positions
                let s1: Vec<_> = above1.iter().map(|&i| rotated1[i]).collect();

                let mut s2: Vec<_> = below2.iter().map(|&i| rotated2[i]).collect();

                s2.extend(s1.iter());

                omol.set_positions(&s2);
                break;
            }
        }

        iloop += 1;
        if iloop >= maxloop {
            bail!("max iterations reached.");
        }
    }

    Ok(omol)
}
// core:1 ends here

// test

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*test][test:1]]
#[test]
fn test_plane_cut_splice() {
    let mut mols1 = gchemol::io::read("tests/files/c6h6-geom1.mol2").expect("geom1");
    let mut mol1 = &mut mols1[0];
    let mut mols2 = gchemol::io::read("tests/files/c6h6-geom2.mol2").expect("geom2");
    let mut mol2 = &mut mols2[0];
    assert_eq!(mol1.natoms(), mol2.natoms());

    let x = plane_cut_and_splice(&mol1, &mol2).expect("plane-cut-and-splice");
    x.to_file("/tmp/aa.xyz").expect("write splice");
}
// test:1 ends here
