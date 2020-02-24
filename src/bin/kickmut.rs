// main

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*main][main:1]]
use gut::prelude::*;
use structopt::*;

use gosh::gchemol;
use gosh::gchemol::prelude::*;
use gosh::gchemol::Molecule;

/// Generate randomly mutated structures
#[derive(Debug, StructOpt)]
struct Cli {
    /// Path to input file containing molecule structure.
    inpfile: String,

    /// Path to output file for writing new molecules.
    outfile: String,

    /// The number of molecules to be generated.
    #[structopt(short = "n", default_value = "1")]
    number_of_molecules: usize,

    /// The number of molecules to be generated.
    #[structopt(short = "m", default_value = "1")]
    mutation_degree: usize,
}

fn main() -> Result<()> {
    let args = Cli::from_args();
    gut::cli::setup_logger();

    let old_mol = Molecule::from_file(args.inpfile)?;
    let mut mols = vec![];
    for _ in 0..args.number_of_molecules {
        let new_mol = kickstart::adhoc::random_bond_mutate(&old_mol, args.mutation_degree).unwrap();
        mols.push(new_mol);
    }

    gchemol::io::write(args.outfile, &mols);
    Ok(())
}
// main:1 ends here
