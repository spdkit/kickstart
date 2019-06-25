// lib.rs
// :PROPERTIES:
// :header-args: :tangle src/lib.rs
// :END:

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*lib.rs][lib.rs:1]]
pub mod kickstart;
pub mod config;
pub mod search;

pub use self::kickstart::kick;

pub(crate) mod common {
    pub use quicli::prelude::*;
    pub type Result<T> = ::std::result::Result<T, Error>;
}
// lib.rs:1 ends here
