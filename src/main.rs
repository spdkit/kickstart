// [[file:../kickstart.note::*imports][imports:1]]
use std::path::PathBuf;
// imports:1 ends here

// [[file:../kickstart.note::ae9bb309][ae9bb309]]
use gut::cli::*;
use gut::prelude::*;

/// Chemical structure explorer
#[derive(Debug, StructOpt)]
#[clap(author, version, about)]
struct Cli {
    /// Prints default configuration.
    #[structopt(long = "print", short = 'p')]
    print: bool,

    /// List calculated items in database.
    #[structopt(long = "list", short = 'l')]
    list: bool,

    /// Run genetic search.
    #[structopt(long = "run", short = 'r')]
    run: bool,
}

fn main() -> Result<()> {
    gut::cli::setup_logger();

    let pkg_version = env!("CARGO_PKG_VERSION");
    let pkg_name = env!("CARGO_PKG_NAME");

    let args = Cli::parse();

    if args.print {
        println!("{:#^72}", " default configuration ");
        kickstart::print_default_config();
    } else if args.list {
        // setup a pager like `less` cmd
        pager::Pager::with_pager("less").setup();
        kickstart::list_db()?;
    } else if args.run {
        let msg = format!(" {} v{} ", pkg_name, pkg_version);
        println!("{:#^80}", msg);
        kickstart::genetic_search()?;
    } else {
        Cli::clap().print_help()?;
    }

    Ok(())
}
// ae9bb309 ends here
