// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use crate::common::*;

use gosh_db::DbConnection;
// imports:1 ends here

// global

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*global][global:1]]
// global database connection
lazy_static! {
    pub(crate) static ref KICKSTART_DB_CONNECTION: DbConnection = {
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
// global:1 ends here
