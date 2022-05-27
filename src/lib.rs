// [[file:../kickstart.note::5bd520e5][5bd520e5]]
#[macro_use]
extern crate lazy_static;

mod config;
mod core;
mod crossover;
mod event_loop;
mod exploitation;
mod exploration;
mod interactive;
mod kickstart;
mod model;
mod mutation;
mod search;

pub mod cli;
// 5bd520e5 ends here

// [[file:../kickstart.note::73ac23dc][73ac23dc]]
pub use gosh::gchemol;

pub use crate::config::print_default_config; // for kickstart/main
pub use crate::core::list_db;
pub use crate::event_loop::enter_main_loop; // for kickstart/main
pub use crate::kickstart::kick; // for bin/kickgen
pub use crate::search::genetic_search; // for kickstart/main

pub(crate) mod common {
    pub use gosh::gchemol;
    pub use gut::prelude::*;
}
// 73ac23dc ends here

// [[file:../kickstart.note::728878b3][728878b3]]
/// Return a version string at compile time
// mhm 0.0.28 (ac3391f 2022-05-26)
pub fn app_version() -> String {
    let version = env!("CARGO_PKG_VERSION");
    let commit_hash = env!("GIT_HASH");
    let commit_date = env!("COMMIT_DATE");
    format!("{version} ({commit_hash} {commit_date})")
}

use gut::prelude::*;
// 728878b3 ends here
