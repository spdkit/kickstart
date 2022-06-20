// [[file:../kickstart.note::957b83c7][957b83c7]]
use gosh::gchemol;
use std::collections::HashMap;

use gchemol::compat::*;
use gchemol::geom::random::{rand_points_within_sphere, rand_rotate};
use gchemol::prelude::*;
use gchemol::{Atom, Molecule};

use crate::common::*;
// 957b83c7 ends here

// [[file:../kickstart.note::6540c145][6540c145]]
/// rotate the molecule in place
fn rotate_molecule(mol: &mut Molecule) -> Result<()> {
    let positions = mol.positions_vec();
    let new = rand_rotate(&positions);
    mol.set_positions(new);

    Ok(())
}

fn combine_fragments_into_one(fragments: &Vec<Molecule>) -> Molecule {
    use educate::prelude::*;

    assert!(!fragments.is_empty(), "empty list of fragments!");

    // collect atoms
    let mut mol = Molecule::new("combined");
    let mut i = 1;
    for f in fragments {
        for (_, a) in f.atoms() {
            mol.add_atom(i, a.clone());
            i += 1;
        }
    }

    trace!("combined {} fragments, {} atoms.", fragments.len(), mol.natoms());

    // FIXME: handle NAN values in position
    let mut nan = false;
    for (_, a) in mol.atoms() {
        let [x, y, z] = a.position();
        if x.is_nan() {
            nan = true;
        }
    }

    mol.educate().context("failed to educate").unwrap();
    mol
}

fn generate_rand_fragments(fragments: &mut Vec<Molecule>, r: f64) -> Result<()> {
    let n = fragments.len();
    let pts = rand_points_within_sphere(r, n);
    trace!("kick {:} fragments, radius = {:.2}", n, r);
    for i in 0..n {
        let mut mol = &mut fragments[i];
        let positions = mol.positions_vec();
        if positions[0][0].is_nan() {
            dbg!(i);
            dbg!(&positions);
        }

        rotate_molecule(&mut mol)?;
        let center = mol.center_of_geometry();
        let mut p = pts[i];
        for j in 0..3 {
            p[j] -= center[j];
        }
        mol.translate(p);
    }

    Ok(())
}

pub(self) fn kickstart(mut mols: &mut Vec<Molecule>, r: f64) -> Result<Vec<Molecule>> {
    generate_rand_fragments(&mut mols, r)?;

    let mol = combine_fragments_into_one(&mols);
    let mols = mol.fragment();

    Ok(mols)
}
// 6540c145 ends here

// [[file:../kickstart.note::03afa91b][03afa91b]]
// FIXME: read formula
pub fn kick(mol: &Molecule) -> Result<Molecule> {
    let mut mols = mol.fragment();
    trace!("kick {} fragments ...", mols.len());
    if mols.len() <= 1 {
        warn!("cannot break molecule into multiple parts!");
        return Ok(mol.to_owned());
    }

    // initial sphere radius
    let mut radius = 8.0;
    loop {
        mols = kickstart(&mut mols, radius)?;
        radius -= 0.5;
        if radius < 1.0 {
            break;
        }
    }

    let mol = combine_fragments_into_one(&mols);

    Ok(mol)
}

/// Create `nbunch` molecules in random configurations using kickstart algorithm in parallel
///
/// # Parameters
///
/// * parent_mol: initial molecule containing multiple fragments (based on connectivity)
/// * nbunch: the number of configurations to be generated
///
pub fn kick_bunch(parent_mol: &Molecule, nbunch: usize) -> Vec<Molecule> {
    debug!("Creating {} molecules using kickstart.", nbunch);
    (0..nbunch)
        .into_par_iter()
        .map(|_| kick(&parent_mol).expect("kick parent_mol"))
        .collect()
}
// 03afa91b ends here

// [[file:../kickstart.note::b63a9460][b63a9460]]
#[test]
fn test_distribute() {
    use gosh::gchemol::prelude::*;

    let filename = "tests/files/c6h6.mol2";
    let mut mol = Molecule::from_file(filename).expect("mol2 file");
    let mol = kick(&mol).expect("kick new mol");
    assert_eq!(12, mol.natoms());
}
// b63a9460 ends here
