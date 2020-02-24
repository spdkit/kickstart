// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use std::path::PathBuf;
// imports:1 ends here

// cmdline

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*cmdline][cmdline:1]]
use gut::prelude::*;
use structopt::*;

/// Chemical structure explorer
#[derive(Debug, StructOpt)]
struct Cli {
    /// Prints default configuration.
    #[structopt(long = "print", short = "p")]
    print: bool,

    /// List calculated items in database.
    #[structopt(long = "list", short = "l")]
    list: bool,

    /// Run genetic search.
    #[structopt(long = "run", short = "r")]
    run: bool,
}

fn main() -> Result<()> {
    gut::cli::setup_logger();

    let pkg_version = env!("CARGO_PKG_VERSION");
    let pkg_name = env!("CARGO_PKG_NAME");

    let args = Cli::from_args();

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
// cmdline:1 ends here
