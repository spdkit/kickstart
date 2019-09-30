// mods

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*mods][mods:1]]
#[macro_use]
extern crate lazy_static;

mod config;
mod kickstart;
mod search;
mod crossover;
mod mutation;
mod model;
// mods:1 ends here

// exports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*exports][exports:1]]
pub use gosh::gchemol;

pub use crate::config::print_default_config;
pub use crate::kickstart::kick;
pub use crate::search::genetic_search;

pub(crate) mod common {
    pub use quicli::prelude::*;
    pub type Result<T> = std::result::Result<T, Error>;
}

pub mod adhoc {
    pub use crate::mutation::*;
}
// exports:1 ends here
