// lib.rs
// :PROPERTIES:
// :header-args: :tangle src/lib.rs
// :END:

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*lib.rs][lib.rs:1]]
pub mod config;
pub mod kickstart;
pub mod search;
pub mod search2;
pub mod annealing;
pub mod reinsert;

pub use crate::kickstart::kick;
pub use crate::search::genetic_search;

pub(crate) mod common {
    pub use quicli::prelude::*;
    pub type Result<T> = std::result::Result<T, Error>;
}
// lib.rs:1 ends here
