// [[file:../kickstart.note::5730952f][5730952f]]
use super::*;

use gosh::gchemol;
use gosh::gchemol::prelude::*;
use gosh::gchemol::Molecule;
use gut::cli::*;
// 5730952f ends here

// [[file:../kickstart.note::2a2bd44d][2a2bd44d]]
/// Generate randomly mutated structures
#[derive(Debug, StructOpt)]
#[clap(author, version, about)]
pub struct KickMut {
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

impl KickMut {
    pub fn enter_main() -> Result<()> {
        let args = Self::parse();
        gut::cli::setup_logger();

        let old_mol = Molecule::from_file(args.inpfile)?;
        let mut mols = vec![];
        for _ in 0..args.number_of_molecules {
            let new_mol = crate::mutation::random_bond_mutate(&old_mol, args.mutation_degree).unwrap();
            mols.push(new_mol);
        }

        gchemol::io::write(args.outfile, &mols);
        Ok(())
    }
}
// 2a2bd44d ends here

// [[file:../kickstart.note::251fd511][251fd511]]
/// Generate random structures
#[derive(Debug, StructOpt)]
#[clap(author, version, about)]
pub struct KickGen {
    /// Path to input file containing molecule structure.
    inpfile: String,

    /// Path to output file for writing new molecules.
    outfile: String,

    /// The number of molecules to be generated.
    #[structopt(short = 'n', default_value = "1")]
    number_of_molecules: usize,
}

impl KickGen {
    pub fn enter_main() -> Result<()> {
        let args = Self::parse();
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
}
// 251fd511 ends here

// [[file:../kickstart.note::ae9bb309][ae9bb309]]
/// Chemical structure explorer
#[derive(Debug, Parser)]
#[clap(author, version, about)]
struct Cli {
    #[clap(flatten)]
    verbose: gut::cli::Verbosity,

    /// Prints default configuration.
    #[clap(long, short)]
    print: bool,

    /// List calculated items in database.
    #[clap(long, short)]
    list: bool,

    /// Sort by energy when listing database.
    #[clap(long)]
    sort: bool,

    /// Run genetic search.
    #[clap(long, short)]
    run: bool,
}

pub fn enter_main() -> Result<()> {
    let args = Cli::parse();
    args.verbose.setup_logger();

    if args.print {
        println!("{:#^72}", " default configuration ");
        crate::print_default_config();
    } else if args.list {
        // setup a pager like `less` cmd
        pager::Pager::with_pager("less").setup();
        crate::list_db(args.sort)?;
    } else if args.run {
        println!("kickstart {:#^80}", app_version());
        crate::genetic_search()?;
    } else {
        Cli::command().print_help()?;
    }

    Ok(())
}
// ae9bb309 ends here
