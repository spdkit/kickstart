// main.rs
// :PROPERTIES:
// :header-args: :tangle src/main.rs
// :END:

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*main.rs][main.rs:1]]
use quicli::prelude::*;
use structopt::*;

use kickstart::config::Config;

/// Chemical structure explorer
#[derive(Debug, StructOpt)]
struct Cli {
    #[structopt(flatten)]
    verbosity: Verbosity,

    /// configuration file for kickstart
    #[structopt(help = "Configuration file in toml format", parse(from_os_str))]
    configfile: std::path::PathBuf,
}

fn main() -> CliResult {
    let args = Cli::from_args();
    args.verbosity.setup_env_logger(&env!("CARGO_PKG_NAME"))?;

    let toml_str = read_file(args.configfile)?;
    let config: Config = toml::from_str(&toml_str)?;

    kickstart::genetic_search(&config)?;

    Ok(())
}
// main.rs:1 ends here
