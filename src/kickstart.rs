// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use std::collections::HashMap;

use gosh::gchemol::geometry::{rand_points_within_sphere, rand_rotate};
use gosh::gchemol::prelude::*;
use gosh::gchemol::{io, Atom, Molecule};

use crate::common::*;
// imports:1 ends here

// new

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*new][new:1]]
/// rotate the molecule in place
fn rotate_molecule(mol: &mut Molecule) -> Result<()> {
    let positions = mol.positions();
    let new = rand_rotate(&positions);
    mol.set_positions(&new)?;

    Ok(())
}

fn combine_fragments_into_one(fragments: &Vec<Molecule>) -> Molecule {
    use educate::prelude::*;

    assert!(!fragments.is_empty(), "empty list of fragments!");

    let mut mol = Molecule::new("combined");
    // collect atoms
    for f in fragments {
        for a in f.atoms() {
            mol.add_atom(a.clone());
        }
    }
    debug!(
        "combined {} fragments, {} atoms.",
        fragments.len(),
        mol.natoms()
    );

    // FIXME: handle NAN values in position
    let mut nan = false;
    for a in mol.atoms() {
        let [x, y, z] = a.position();
        if x.is_nan() {
            nan = true;
        }
    }

    mol.educate().with_context(|_| {
        format!("failed to educate")
    }).unwrap();

    mol
}

fn generate_rand_fragments(fragments: &mut Vec<Molecule>, r: f64) -> Result<()> {
    let n = fragments.len();
    let pts = rand_points_within_sphere(r, n);
    debug!("kick {:} fragments, radius = {:.2}", n, r);
    for i in 0..n {
        let mut mol = &mut fragments[i];
        let positions = mol.positions();
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

fn kickstart(mut mols: &mut Vec<Molecule>, r: f64) -> Result<Vec<Molecule>> {
    generate_rand_fragments(&mut mols, r)?;

    let mol = combine_fragments_into_one(&mols);
    let mols = mol.fragment();

    Ok(mols)
}
// new:1 ends here

// pub

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*pub][pub:1]]
// FIXME: read formula
pub fn kick(mol: &Molecule) -> Result<Molecule> {
    let mut mols = mol.fragment();
    debug!("kick {} fragments ...", mols.len());
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
// pub:1 ends here

// test

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*test][test:1]]
#[test]
fn test_distribute() {
    use gosh::gchemol::prelude::*;

    let filename = "tests/files/c6h6.mol2";
    let mut mol = Molecule::from_file(filename).expect("mol2 file");
    let mol = kick(&mol).expect("kick new mol");
    assert_eq!(12, mol.natoms());
    mol.to_file("/tmp/k.mol2");
}
// test:1 ends here
