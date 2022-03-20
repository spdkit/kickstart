// [[file:../../kickstart.note::251fd511][251fd511]]
use gut::cli::*;
use gut::prelude::*;

use gosh::gchemol;
use gosh::gchemol::prelude::*;
use gosh::gchemol::Molecule;

/// Generate random structures
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
}

fn main() -> Result<()> {
    let args = Cli::parse();
    gut::cli::setup_logger();

    let mol = Molecule::from_file(args.inpfile)?;
    let mut mols = vec![];
    for _ in 0..args.number_of_molecules {
        let mol = kickstart::kick(&mol)?;
        mols.push(mol);
    }

    gchemol::io::write(args.outfile, &mols);
    Ok(())
}
// 251fd511 ends here
