// mods

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*mods][mods:1]]
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
// mods:1 ends here

// exports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*exports][exports:1]]
pub use gosh::gchemol;

pub use crate::config::print_default_config; // for kickstart/main
pub use crate::core::list_db;
pub use crate::event_loop::enter_main_loop; // for kickstart/main
pub use crate::kickstart::kick; // for bin/kickgen
pub use crate::search::genetic_search; // for kickstart/main // for main

pub(crate) mod common {
    pub use quicli::prelude::*;
    pub type Result<T> = std::result::Result<T, Error>;
    pub use gosh::gchemol;
}

pub mod adhoc {
    pub use crate::mutation::*;
}
// exports:1 ends here
