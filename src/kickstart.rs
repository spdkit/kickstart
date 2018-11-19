// new

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*new][new:1]]
use std::collections::HashMap;
use quicli::prelude::*;

use gchemol::{
    Atom,
    Molecule,
    io,
};

use gchemol::geometry::{
    rand_rotate,
    rand_points_within_sphere,
};

/// rotate the molecule in place
fn rotate_molecule(mol: &mut Molecule) {
    let positions = mol.positions();
    let new = rand_rotate(&positions);
    mol.set_positions(&new).expect("assign new positions");
}

// fn create_c6h6() -> Vec<Molecule> {
//     let mut mols = vec![];
//     for e in vec!["C", "H"] {
//         for i in 0..6 {
//             let mut mol = Molecule::new(&format!("{}{}", e, i));
//             mol.add_atom(Atom::new(e, [0.0; 3]));
//             mols.push(mol);
//         }
//     }

//     mols
// }

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
    mol.clean();
    mol
}

fn generate_rand_fragments(fragments: &mut Vec<Molecule>, r: f64) {
    let n = fragments.len();
    let pts = rand_points_within_sphere(r, n);
    println!("kick {:} fragments, radius = {:.2}", n, r);
    for i in 0..n {
        let mut mol = &mut fragments[i];
        rotate_molecule(&mut mol);
        let center = mol.center_of_geometry();
        let mut p = pts[i];
        for j in 0..3 {
            p[j] -= center[j];
        }
        mol.translate(p);
    }
}

fn kickstart(mut mols: &mut Vec<Molecule>, r: f64) -> Vec<Molecule>{
    {
        generate_rand_fragments(&mut mols, r);
    }
    let mol = combine_fragments_into_one(&mols);
    // let filename = format!("/tmp/{:}.mol2", r.round());
    // mol.to_file(filename);
    mol.fragment()
}

// FIXME: read formula
pub fn kick(mol: &Molecule) -> Result<Molecule> {
    info!("kick molecule...");
    let mut mols = mol.fragment();
    if mols.len() <= 1 {
        bail!("cannot break molecule into multiple parts!");
    }

    // initial sphere radius
    let mut radius = 5.0;
    loop {
        mols = kickstart(&mut mols, radius);
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
    use gchemol::prelude::*;

    let filename = "tests/files/c6h6.mol2";
    let mut mol = Molecule::from_file(filename).expect("mol2 file");
    let mol = kick(&mol).expect("kick new mol");
    assert_eq!(12, mol.natoms());
    mol.to_file("/tmp/k.mol2");
}
// new:1 ends here
