// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use genevo::prelude::*;

use gosh::gchemol::prelude::*;
use gosh::gchemol::{io, Atom, Molecule};

use crate::common::*;
use crate::config::Config;
// imports:1 ends here

// base

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*base][base:1]]
/// The Genotype
#[derive(Clone, Debug, PartialEq)]
pub struct MoleculeGenome {
    name: String,
    data: Vec<(usize, [f64; 3])>,
}

impl Genotype for MoleculeGenome {
    type Dna = Vec<(usize, [f64; 3])>;
}

pub trait ToGenome {
    fn to_genome(&self) -> MoleculeGenome;
}

/// Create a random string of length `n` for naming a genome
fn random_name(n: usize) -> String {
    use rand::distributions::Alphanumeric;
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    rng.sample_iter(&Alphanumeric).take(n).collect()
}

impl ToGenome for Molecule {
    fn to_genome(&self) -> MoleculeGenome {
        let mut g = vec![];
        for a in self.sorted().atoms() {
            let n = a.number();
            let p = a.position();
            g.push((n, p));
        }

        MoleculeGenome {
            name: random_name(5),
            data: g,
        }
    }
}

impl MoleculeGenome {
    pub fn to_molecule(&self) -> Molecule {
        let mut mol = Molecule::new(&self.name);

        for (n, [x, y, z]) in &self.data {
            let a = Atom::build().element(*n).position(*x, *y, *z).finish();
            mol.add_atom(a);
        }

        mol.rebond();
        mol
    }
}
// base:1 ends here
