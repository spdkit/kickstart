// [[file:../kickstart.note::5712699d][5712699d]]
use crate::common::*;
use gosh::gchemol;

use gchemol::compat::*;
use gchemol::{Atom, Molecule};
use spdkit::prelude::*;
// 5712699d ends here

// [[file:../kickstart.note::62013f9d][62013f9d]]
use vecfx::*;

type OF64 = vecfx::OrderedFloat<f64>;

/// The Genotype for molecule. Core data structure for evolution.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MolGenome {
    name: String,
    age: usize,
    data: Vec<(usize, [OF64; 3])>,
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
    /// Uniq genome ID
    pub fn uid(&self) -> &str {
        &self.name
    }

    /// Healthy point.
    pub fn hp(&self) -> f64 {
        (-BETA_FACTOR * self.age as f64).exp()
    }

    /// Return a cloned MolGenome with larger age.
    pub fn aged(&self) -> Self {
        let mut new = self.clone();
        new.age += 1;
        new
    }
}
// 62013f9d ends here

// [[file:../kickstart.note::807e3191][807e3191]]
use std::hash::{Hash, Hasher};

impl Hash for MolGenome {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // self.name.hash(state);
        self.data.hash(state);
    }
}

impl PartialEq for MolGenome {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl Eq for MolGenome {}
// 807e3191 ends here

// [[file:../kickstart.note::47e1ae28][47e1ae28]]
/// The evaluated energy with molecule structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EvaluatedGenome {
    pub genome: MolGenome,
    pub energy: f64,
}

impl EvaluatedGenome {
    /// unique ID for saving into and retrieving from database.
    fn uid(&self) -> &str {
        self.genome.uid()
    }
}
// 47e1ae28 ends here

// [[file:../kickstart.note::376404cc][376404cc]]
use self::db::KICKSTART_DB_CONNECTION as Db;
use gosh::db::prelude::*;

impl EvaluatedGenome {
    /// Return total number of evaluations
    pub fn number_of_evaluations() -> usize {
        Self::collection_size(&Db).expect("db: list failure") as usize
    }

    /// List items found in checkpointing database
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
    pub fn energy(&self) -> f64 {
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

// [[file:../kickstart.note::e4f2cc3b][e4f2cc3b]]
const GENOME_NAME_LENGTH: usize = 8;

pub trait ToGenome {
    fn encode(&self) -> MolGenome;

    // FIXME: adhoc
    fn encode_as_evaluated(&self) -> EvaluatedGenome {
        unimplemented!()
    }
}

impl MolGenome {
    /// Re-create molecule from MolGenome
    pub fn decode(&self) -> Molecule {
        use educate::prelude::*;

        let mut mol = Molecule::from_atoms(self.data.iter().map(|&(sym, coords)| (sym, coords.map(|x| x.into()))));

        // recreate connectivity from current positions
        mol.educated_rebond().unwrap();
        mol.set_title(&self.name);
        mol
    }
}
// e4f2cc3b ends here

// [[file:../kickstart.note::6c1c6f41][6c1c6f41]]
/// Create a random string of length `n` for naming a genome
fn random_name(n: usize) -> String {
    use rand::distributions::Alphanumeric;

    let mut rng = thread_rng();
    rng.sample_iter(&Alphanumeric).take(n).map(char::from).collect()
}

mod hash {
    use super::*;
    use sha2::{Digest, Sha256};

    fn create_hash(msg: &str) -> String {
        let mut hasher = Sha256::new();
        hasher.update(msg);
        format!("{:x}", hasher.finalize()).chars().take(8).collect()
    }

    impl MolGenome {
        pub fn encode_molecule(mol: &Molecule) -> Self {
            let mut g = vec![];
            for (_, a) in mol.sorted().atoms() {
                let n = a.number();
                let p = a.position().map(|x| x.as_ordered_float());
                g.push((n, p));
            }

            let mut mol = mol.clone();
            mol.rebond();
            let fp = mol.fingerprint();
            let name = self::create_hash(&fp);

            Self { name, age: 0, data: g }
        }
    }

    #[test]
    #[ignore]
    fn test_hash_code() {
        dbg!(create_hash("xxit"));
    }
}
// 6c1c6f41 ends here

// [[file:../kickstart.note::f51a8eee][f51a8eee]]
use std::sync::atomic;

pub type JobFlag = atomic::AtomicUsize;

#[derive(Eq, PartialEq, Debug)]
pub enum JobType {
    /// Run program
    Run,
    /// Stop program.
    Stop,
    /// Edit internal state.
    Edit,
}

impl JobType {
    pub fn from(flag: &JobFlag) -> Self {
        match flag.load(atomic::Ordering::SeqCst) {
            0 => Self::Run,
            1 => Self::Stop,
            2 => Self::Edit,
            _ => unreachable!(),
        }
    }

    pub fn flag(&self) -> usize {
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
// f51a8eee ends here

// [[file:../kickstart.note::44895e0d][44895e0d]]
use gosh::model::*;

/// avoid recalculation with database caching
#[derive(Debug, Clone)]
pub struct MolIndividual;

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
            genome: MolGenome::encode_molecule(mol),
            energy,
        };
        let key = evaluated.uid();
        trace!("saving result with key {}", key);
        evaluated.put_into_collection(&Db, key).expect("db write failure");

        evaluated
    }
}

/// Return the number of evaluation of molecules
pub fn get_number_of_evaluations() -> usize {
    EvaluatedGenome::number_of_evaluations()
}
// 44895e0d ends here
