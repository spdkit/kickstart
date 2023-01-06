// [[file:../kickstart.note::5bd520e5][5bd520e5]]
#[macro_use]
extern crate lazy_static;

mod config;
mod core;
mod crossover;
// mod event_loop;
mod exploitation;
mod exploration;
// mod interactive;
mod kickstart;
mod model;
mod mutation;
mod search;

mod hardsphere;

pub mod cli;
// 5bd520e5 ends here

// [[file:../kickstart.note::73ac23dc][73ac23dc]]
use gosh::gchemol;
use gut::prelude::*;

use crate::config::print_default_config; // for kickstart/main
use crate::core::list_db;
use crate::search::genetic_search; // for kickstart/main

mod common {
    pub use gosh::gchemol;
    pub use gut::prelude::*;

    pub use std::path::{Path, PathBuf};
}
// 73ac23dc ends here

// [[file:../kickstart.note::728878b3][728878b3]]
/// Return a version string at compile time
// mhm 0.0.28 (ac3391f 2022-05-26)
fn app_version() -> String {
    let version = env!("CARGO_PKG_VERSION");
    let commit_hash = env!("GIT_HASH");
    let commit_date = env!("COMMIT_DATE");
    format!("{version} ({commit_hash} {commit_date})")
}
// 728878b3 ends here

// [[file:../kickstart.note::58ff7b77][58ff7b77]]
#[cfg(feature = "adhoc")]
/// Docs for local mods
pub mod docs {
    macro_rules! export_doc {
        ($l:ident) => {
            pub mod $l {
                pub use crate::$l::*;
            }
        };
    }

    export_doc!(core);
    export_doc!(exploitation);
    export_doc!(config);
    export_doc!(hardsphere);
}
// 58ff7b77 ends here
