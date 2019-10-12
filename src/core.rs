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

impl MolGenome {
    pub(crate) fn uid(&self) -> &str {
        &self.name
    }
}
// genome:1 ends here

// evaluated

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*evaluated][evaluated:1]]
/// The evaluated energy with molecule structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct EvaluatedGenome {
    pub genome: MolGenome,
    pub energy: f64,
}

impl EvaluatedGenome {
    /// unique ID for saving into and retrieving from database.
    fn uid(&self) -> &str {
        self.genome.uid()
    }
}
// evaluated:1 ends here

// database

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*database][database:1]]
use self::db::KICKSTART_DB_CONNECTION as Db;
use gosh_db::prelude::*;

impl Collection for EvaluatedGenome {
    fn collection_name() -> String {
        "EvaluatedGenome".into()
    }
}

impl EvaluatedGenome {
    pub(crate) fn number_of_evaluations() -> usize {
        Self::collection_size(&Db).expect("db: list failure") as usize
    }
}

impl MolGenome {
    /// Retrieve energy from db
    pub(crate) fn energy(&self) -> f64 {
        let key = self.uid();
        let evaluated =
            EvaluatedGenome::get_from_collection(&Db, key).expect("db: read energy failure");
        evaluated.energy
    }
}

// global database connection
mod db {
    use crate::common::*;
    use gosh_db::prelude::*;
    use gosh_db::DbConnection;

    lazy_static! {
        pub(super) static ref KICKSTART_DB_CONNECTION: DbConnection = {
            let dbvar = "GOSH_DATABASE_URL";
            let default_db = format!("{}.db", env!("CARGO_PKG_NAME"));
            if std::env::var(dbvar).is_err() {
                info!("Use default db file: {}", default_db);
                std::env::set_var(dbvar, default_db);
            }
            let db = DbConnection::establish().expect("gosh db");
            db
        };
    }
}
// database:1 ends here

// genome/molecule mapping
// genotype <=> phenotype conversion

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*genome/molecule%20mapping][genome/molecule mapping:1]]
const GENOME_NAME_LENGTH: usize = 8;

pub(crate) trait ToGenome {
    fn encode(&self) -> MolGenome;

    // FIXME: adhoc
    fn encode_as_evaluated(&self) -> EvaluatedGenome {
        unimplemented!()
    }
}

impl MolGenome {
    /// Re-create molecule from MolGenome
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
    let mut rng = thread_rng();
    rng.sample_iter(&Alphanumeric).take(n).collect()
}

fn encode_molecule(mol: &Molecule) -> MolGenome {
    let mut g = vec![];
    for a in mol.sorted().atoms() {
        let n = a.number();
        let p = a.position();
        g.push((n, p));
    }

    MolGenome {
        name: random_name(GENOME_NAME_LENGTH),
        data: g,
    }
}
// genome/molecule mapping:1 ends here

// public

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*public][public:1]]
use gosh::models::*;

/// avoid recalculation with database caching
#[derive(Debug, Clone)]
pub(crate) struct MolIndividual;

/// for Valuer: valuer.create_individuals
impl EvaluateObjectiveValue<MolGenome> for MolIndividual {
    fn evaluate(&self, genome: &MolGenome) -> f64 {
        genome.energy()
    }
}

/// Encode computed molecule as `MolGenome` for evolution. The computed
/// results will be cached in database for spdkit later retrieving.
impl ToGenome for ModelProperties {
    fn encode(&self) -> MolGenome {
        self.encode_as_evaluated().genome
    }

    fn encode_as_evaluated(&self) -> EvaluatedGenome {
        let energy = self.energy.expect("no energy");
        let mol = self.molecule.as_ref().expect("no molecule");
        let evaluated = EvaluatedGenome {
            genome: encode_molecule(mol),
            energy,
        };
        let key = evaluated.uid();
        trace!("saving result with key {}", key);
        evaluated
            .put_into_collection(&Db, key)
            .expect("db write failure");

        evaluated
    }
}

/// Return the number of evaluation of molecules
pub(crate) fn get_number_of_evaluations() -> usize {
    EvaluatedGenome::number_of_evaluations()
}
// public:1 ends here
