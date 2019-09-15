// lib.rs
// :PROPERTIES:
// :header-args: :tangle src/lib.rs
// :END:

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*lib.rs][lib.rs:1]]
#[macro_use]
extern crate lazy_static;

mod search;
mod config;
mod kickstart;
mod annealing;

pub use crate::kickstart::kick;
pub use crate::search::genetic_search;
pub use crate::config::print_default_config;

pub(crate) mod common {
    pub use quicli::prelude::*;
    pub type Result<T> = std::result::Result<T, Error>;
}
// lib.rs:1 ends here
