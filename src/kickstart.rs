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
    let mut mol = Molecule::new("combined");
    // collect atoms
    for f in fragments {
        for a in f.atoms() {
            mol.add_atom(a.clone());
        }
    }
    println!("combined {} fragments.", fragments.len());

    // perceive bonding connectivity
    mol.rebond();

    let mut nan = false;
    for a in mol.atoms() {
        let [x, y, z] = a.position();
        if x.is_nan() {
            nan = true;
        }
    }

    // if !nan {
    //     mol.to_file("/tmp/nan1.mol2");
    // }

    mol.clean().expect("clean");
    mol
}

fn generate_rand_fragments(fragments: &mut Vec<Molecule>, r: f64) -> Result<()> {
    let n = fragments.len();
    let pts = rand_points_within_sphere(r, n);
    println!("kick {:} fragments, radius = {:.2}", n, r);
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

// FIXME: read formula
pub fn kick(mol: &Molecule) -> Result<Molecule> {
    let mut mols = mol.fragment();
    info!("kick {} fragments ...", mols.len());
    if mols.len() <= 1 {
        bail!("cannot break molecule into multiple parts!");
    }

    // initial sphere radius
    let mut radius = 6.0;
    loop {
        mols = kickstart(&mut mols, radius)?;
        radius -= 0.5;
        if radius < 1.5 {
            break;
        }
    }

    let mol = combine_fragments_into_one(&mols);

    Ok(mol)
}

#[test]
fn test_distribute() {
    use gosh::gchemol::prelude::*;

    let filename = "tests/files/c6h6.mol2";
    let mut mol = Molecule::from_file(filename).expect("mol2 file");
    let mol = kick(&mol).expect("kick new mol");
    assert_eq!(12, mol.natoms());
    mol.to_file("/tmp/k.mol2");
}
// new:1 ends here
