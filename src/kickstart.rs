// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::86373bab-a76e-4baf-842b-dc6cebc0bed7][86373bab-a76e-4baf-842b-dc6cebc0bed7]]
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
    mol.set_positions(new).expect("assign new positions");
}

fn create_c6h6() -> Vec<Molecule> {
    let mut mols = vec![];
    for e in vec!["C", "H"] {
        for i in 0..6 {
            let mut mol = Molecule::new(&format!("{}{}", e, i));
            mol.add_atom(Atom::new(e, [0.0; 3]));
            mols.push(mol);
        }
    }

    mols
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
    let filename = format!("/tmp/{:}.mol2", r.round());
    mol.to_file(filename);

    mol.fragment()
}

#[test]
fn test_distribute() {
    let mut mols = create_c6h6();
    let mut mols = kickstart(&mut mols, 5.0);
    let mut mols = kickstart(&mut mols, 4.0);
    let mut mols = kickstart(&mut mols, 3.0);
    let mut mols = kickstart(&mut mols, 2.0);
    let mut mols = kickstart(&mut mols, 1.0);
}
// 86373bab-a76e-4baf-842b-dc6cebc0bed7 ends here
