// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use crate::common::*;

use gchemol::{Atom, Molecule};
use gosh::gchemol;
use spdkit::prelude::*;
// imports:1 ends here

// genome

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*genome][genome:1]]
/// The Genotype for molecule
#[derive(Clone, Debug, Serialize, Deserialize)]
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

// evaluated genome

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*evaluated%20genome][evaluated genome:1]]
/// The evaluated energy with molecule structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct EvaluatedGenome {
    genome: MolGenome,
    energy: f64,
}

impl EvaluatedGenome {
    /// unique ID for saving into and retrieving from database.
    fn uid(&self) -> &str {
        &self.genome.name
    }
}

/// Convenient methods for accessing attributes
impl EvaluatedGenome {
    pub(crate) fn energy(&self) -> f64 {
        self.energy
    }

    pub(crate) fn genome(&self) -> &MolGenome {
        &self.genome
    }
}
// evaluated genome:1 ends here

// base

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*base][base:1]]
use crate::database::KICKSTART_DB_CONNECTION as Db;

use gosh_db::prelude::*;

impl Collection for EvaluatedGenome {
    fn collection_name() -> String {
        "EvaluatedGenome".into()
    }
}
// base:1 ends here

// evaluation

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*evaluation][evaluation:1]]
// impl EvaluatedGenome {
//     /// Save calculated result into default database.
//     pub(crate) fn save_into_database(&self, genome: &MolGenome, energy: f64) -> Result<()> {
//         let key = &genome.name;
//         trace!("saving result with key {}", key);
//         self.put_into_collection(&Db, key)?;

//         Ok(())
//     }
// }

impl EvaluateObjectiveValue<MolGenome> for CachedMolIndividual {
    fn evaluate(&self, genome: &MolGenome) -> f64 {
        let evaluated = self.evaluate_genome(genome).unwrap();
        evaluated.energy
    }
}

/// avoid recalculation with database caching
#[derive(Debug, Clone)]
pub(crate) struct CachedMolIndividual;

impl CachedMolIndividual {
    /// Evaluate with caching using database.
    fn evaluate_genome(&self, genome: &MolGenome) -> Result<EvaluatedGenome> {
        let key = &genome.name;
        match EvaluatedGenome::get_from_collection(&Db, key) {
            Ok(evaluated) => Ok(evaluated),
            // FIXME: handle not-found error
            Err(e) => {
                let evaluated = self.evaluate_new(genome)?;
                evaluated.put_into_collection(&Db, key)?;
                Ok(evaluated)
            }
        }
    }

    /// Evaluate new structure.
    fn evaluate_new(&self, genome: &MolGenome) -> Result<EvaluatedGenome> {
        use crate::model::*;

        let energy = if let Ok(energy) = genome.decode().get_energy() {
            info!("evaluated indv {}, energy = {:-12.5}", genome, energy);
            energy
        } else {
            warn!("Calculation failure. No energy found for {}", genome);
            std::f64::MAX
        };

        let evaluated = EvaluatedGenome {
            genome: genome.clone(),
            energy,
        };
        Ok(evaluated)
    }
}
// evaluation:1 ends here

// genome/molecule mapping
// genotype <=> phenotype conversion

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*genome/molecule%20mapping][genome/molecule mapping:1]]
const GENOME_NAME_LENGTH: usize = 8;

pub(crate) trait ToGenome {
    fn encode(&self) -> MolGenome;
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
            name: random_name(GENOME_NAME_LENGTH),
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

/// Create a random string of length `n` for naming a genome
fn random_name(n: usize) -> String {
    use rand::distributions::Alphanumeric;
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    rng.sample_iter(&Alphanumeric).take(n).collect()
}
// genome/molecule mapping:1 ends here
