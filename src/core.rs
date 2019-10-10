// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use crate::common::*;

use gchemol::{Atom, Molecule};
use gosh::gchemol;
use spdkit::prelude::*;
// imports:1 ends here

// genome

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*genome][genome:1]]
#[derive(Debug, Clone)]
pub(crate) struct MolIndividual;

/// The Genotype for molecule
#[derive(Clone, Debug)]
pub(crate) struct MolGenome {
    name: String,
    data: Vec<(usize, [f64; 3])>,
}

impl std::fmt::Display for MolGenome {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:}", self.name)
    }
}

impl spdkit::individual::Genome for MolGenome {}
// genome:1 ends here

// genome/molecule mapping
// genotype <=> phenotype conversion

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*genome/molecule%20mapping][genome/molecule mapping:1]]
pub(crate) trait ToGenome {
    fn encode(&self) -> MolGenome;
}

/// Create a random string of length `n` for naming a genome
fn random_name(n: usize) -> String {
    use rand::distributions::Alphanumeric;
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    rng.sample_iter(&Alphanumeric).take(n).collect()
}

impl ToGenome for Molecule {
    fn encode(&self) -> MolGenome {
        let mut g = vec![];
        for a in self.sorted().atoms() {
            let n = a.number();
            let p = a.position();
            g.push((n, p));
        }

        MolGenome {
            name: random_name(5),
            data: g,
        }
    }
}

impl MolGenome {
    pub(crate) fn decode(&self) -> Molecule {
        use educate::prelude::*;

        let mut mol = Molecule::new(&self.name);

        for (n, [x, y, z]) in &self.data {
            let a = Atom::build().element(*n).position(*x, *y, *z).finish();
            mol.add_atom(a);
        }

        mol.educated_rebond().unwrap();
        mol
    }
}
// genome/molecule mapping:1 ends here
