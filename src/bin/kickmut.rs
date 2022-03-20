// [[file:../../kickstart.note::b048123f][b048123f]]
use gut::cli::*;
use gut::prelude::*;

use gosh::gchemol;
use gosh::gchemol::prelude::*;
use gosh::gchemol::Molecule;

/// Generate randomly mutated structures
#[derive(Debug, StructOpt)]
#[clap(author, version, about)]
struct Cli {
    /// Path to input file containing molecule structure.
    inpfile: String,

    /// Path to output file for writing new molecules.
    outfile: String,

    /// The number of molecules to be generated.
    #[structopt(short = 'n', default_value = "1")]
    number_of_molecules: usize,

    /// The number of molecules to be generated.
    #[structopt(short = 'm', default_value = "1")]
    mutation_degree: usize,
}

fn main() -> Result<()> {
    let args = Cli::parse();
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
// b048123f ends here
