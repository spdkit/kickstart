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
const BETA_FACTOR: f64 = 0.1;

/// The Genotype for molecule. Core data structure for evolution.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MolGenome {
    fp: String,
    age: usize,
    data: Vec<(usize, [OF64; 3])>,
}

impl std::fmt::Display for MolGenome {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        // output name and current healthy point
        write!(f, "{:} (hp = {:.2})", self.uid(), self.hp())
    }
}

impl spdkit::individual::Genome for MolGenome {}

impl MolGenome {
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
        self.data.hash(state);
    }
}

impl PartialEq for MolGenome {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl Eq for MolGenome {}

impl MolGenome {
    /// Uniq genome ID
    pub fn uid(&self) -> &str {
        // gut::utils::hash_code(&self.data)
        &self.fp
    }

    /// Enode `Molecule` into core data structure for evolution
    pub fn encode_from_molecule(mol: &Molecule) -> Self {
        let mut g = vec![];
        let mut mol = mol.clone();
        mol.rebond();
        // mol.reorder_cannonically();
        for (_, a) in mol.atoms() {
            let n = a.number();
            let p = a.position().map(|x| x.as_ordered_float());
            g.push((n, p));
        }
        let fp = mol.fingerprint();

        Self { age: 0, data: g, fp }
    }
}
// 807e3191 ends here

// [[file:../kickstart.note::376404cc][376404cc]]
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

/// list calculated data in checkpoint database for cli uses
pub fn list_db(sort: bool) -> Result<()> {
    EvaluatedGenome::list_db(sort)
}

// global database connection
mod db {
    use super::*;
    use gosh::db::prelude::*;
    use gosh::db::DbConnection;

    lazy_static! {
        static ref KICKSTART_DB_CONNECTION: DbConnection = {
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

    use KICKSTART_DB_CONNECTION as Db;
    impl EvaluatedGenome {
        /// Return total number of evaluations
        pub fn number_of_evaluations() -> usize {
            Self::collection_size(&Db).expect("db: list failure") as usize
        }

        /// Put evaluated data into database
        pub(super) fn put_into_db(&self) -> Result<()> {
            let key = self.uid();
            trace!("saving result with key {}", key);
            // it happens, especially when run in parallel
            // it is safe to ignore
            if let Err(err) = self.put_into_collection(&Db, key) {
                warn!("write db failure: {err:?}");
            }
            Ok(())
        }

        /// List items found in checkpointing database
        ///
        /// # Parameters
        ///
        /// * sort: sort items by energy
        pub fn list_db(sort: bool) -> Result<()> {
            let mut items = Self::list_collection(&Db)?;
            if items.is_empty() {
                error!("No items in db.");
            } else {
                println!("Found {} items.", items.len());
                let width = items[0].uid().len();
                println!("{:^width$} => {:^12}", "key", "energy");

                if sort {
                    items.sort_by_key(|a| a.energy.as_ordered_float());
                }
                for eg in items {
                    let key = eg.uid();
                    println!("{:^width$} => {:<-12.4}", key, eg.energy);
                }
            }
            Ok(())
        }
    }

    impl MolGenome {
        /// Retrieve energy from db
        pub fn get_energy(&self) -> Option<f64> {
            let key = self.uid();
            EvaluatedGenome::get_from_collection(&Db, &key).ok().map(|item| item.energy)
        }
    }
}
// 376404cc ends here

// [[file:../kickstart.note::e4f2cc3b][e4f2cc3b]]
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
        mol.set_title(self.uid());
        mol
    }
}
// e4f2cc3b ends here

// [[file:../kickstart.note::44895e0d][44895e0d]]
use gosh::model::*;

/// avoid recalculation with database caching
#[derive(Debug, Clone)]
pub struct MolIndividual;

/// for Valuer: valuer.create_individuals
impl EvaluateObjectiveValue<MolGenome> for MolIndividual {
    fn evaluate(&self, genome: &MolGenome) -> f64 {
        // FIXME: review required
        genome.get_energy().unwrap()
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
            genome: MolGenome::encode_from_molecule(mol),
            energy,
        };
        evaluated.put_into_db();

        evaluated
    }
}

/// Return the number of evaluation of molecules
pub fn get_number_of_evaluations() -> usize {
    EvaluatedGenome::number_of_evaluations()
}
// 44895e0d ends here
