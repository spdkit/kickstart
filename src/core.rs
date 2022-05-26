// [[file:../kickstart.note::*imports][imports:1]]
use crate::common::*;
use gosh::gchemol;

use gchemol::compat::*;
use gchemol::{Atom, Molecule};
use spdkit::prelude::*;
// imports:1 ends here

// [[file:../kickstart.note::*genome][genome:1]]
/// The Genotype for molecule
#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct MolGenome {
    name: String,
    age: usize,
    data: Vec<(usize, [f64; 3])>,
}

impl std::fmt::Display for MolGenome {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        // output name and current healthy point
        write!(f, "{:} (hp = {:.2})", self.name, self.hp())
    }
}

impl spdkit::individual::Genome for MolGenome {}

const BETA_FACTOR: f64 = 0.1;
impl MolGenome {
    pub(crate) fn uid(&self) -> &str {
        &self.name
    }

    /// Healthy point.
    pub(crate) fn hp(&self) -> f64 {
        (-BETA_FACTOR * self.age as f64).exp()
    }

    /// Return a cloned MolGenome with larger age.
    pub(crate) fn aged(&self) -> Self {
        let mut new = self.clone();
        new.age += 1;
        new
    }
}
// genome:1 ends here

// [[file:../kickstart.note::*evaluated][evaluated:1]]
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

// [[file:../kickstart.note::376404cc][376404cc]]
use self::db::KICKSTART_DB_CONNECTION as Db;
use gosh::db::prelude::*;

// impl Collection for EvaluatedGenome {
//     fn collection_name() -> String {
//         "EvaluatedGenome".into()
//     }
// }

impl EvaluatedGenome {
    pub(crate) fn number_of_evaluations() -> usize {
        Self::collection_size(&Db).expect("db: list failure") as usize
    }

    pub fn list_db() -> Result<()> {
        let mut items = Self::list_collection(&Db)?;
        if items.is_empty() {
            error!("No items in db.");
        } else {
            println!("Found {} items.", items.len());
            println!("{:^width$} => {:^12}", "key", "energy", width = items[0].uid().len());

            items.sort_by(|a, b| a.energy.partial_cmp(&b.energy).unwrap_or(std::cmp::Ordering::Less));
            for eg in items {
                let key = eg.uid();
                println!("{} => {:<-12.4}", key, eg.energy);
            }
        }
        Ok(())
    }
}

pub fn list_db() -> Result<()> {
    EvaluatedGenome::list_db()
}

impl MolGenome {
    /// Retrieve energy from db
    pub(crate) fn energy(&self) -> f64 {
        let key = self.uid();
        let evaluated = EvaluatedGenome::get_from_collection(&Db, key).expect("db: read energy failure");
        evaluated.energy
    }
}

// global database connection
mod db {
    use crate::common::*;
    use gosh::db::prelude::*;
    use gosh::db::DbConnection;

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
// 376404cc ends here

// [[file:../kickstart.note::*genome/molecule mapping][genome/molecule mapping:1]]
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

        let mut mol = Molecule::from_atoms(self.data.clone());
        mol.set_title(&self.name);

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
    for (_, a) in mol.sorted().atoms() {
        let n = a.number();
        let p = a.position();
        g.push((n, p));
    }

    MolGenome {
        name: random_name(GENOME_NAME_LENGTH),
        age: 0,
        data: g,
    }
}
// genome/molecule mapping:1 ends here

// [[file:../kickstart.note::*job control][job control:1]]
use std::sync::atomic;

pub(crate) type JobFlag = atomic::AtomicUsize;

#[derive(Eq, PartialEq, Debug)]
pub(crate) enum JobType {
    /// Run program
    Run,
    /// Stop program.
    Stop,
    /// Edit internal state.
    Edit,
}

impl JobType {
    pub(crate) fn from(flag: &JobFlag) -> Self {
        match flag.load(atomic::Ordering::SeqCst) {
            0 => Self::Run,
            1 => Self::Stop,
            2 => Self::Edit,
            _ => unreachable!(),
        }
    }

    pub(crate) fn flag(&self) -> usize {
        match self {
            Self::Run => 0,
            Self::Stop => 1,
            Self::Edit => 2,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::sync::atomic;

    #[test]
    fn test_job_flag() {
        let flag = JobFlag::new(JobType::Run.flag());
        dbg!(flag);
    }
}
// job control:1 ends here

// [[file:../kickstart.note::*public][public:1]]
use gosh::model::*;

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
        let energy = self.get_energy().expect("no energy");
        let mol = self.get_molecule().as_ref().expect("no molecule").clone();
        let evaluated = EvaluatedGenome {
            genome: encode_molecule(mol),
            energy,
        };
        let key = evaluated.uid();
        trace!("saving result with key {}", key);
        evaluated.put_into_collection(&Db, key).expect("db write failure");

        evaluated
    }
}

/// Return the number of evaluation of molecules
pub(crate) fn get_number_of_evaluations() -> usize {
    EvaluatedGenome::number_of_evaluations()
}
// public:1 ends here
