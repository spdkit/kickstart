// lib.rs
// :PROPERTIES:
// :header-args: :tangle src/lib.rs
// :END:


// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*lib.rs][lib.rs:1]]
#[macro_use] extern crate quicli;
#[macro_use] extern crate approx;
extern crate rand;
extern crate itertools;
extern crate petgraph;

extern crate gchemol;

pub mod kickstart;
pub use self::kickstart::kick;
// lib.rs:1 ends here
