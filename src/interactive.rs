// [[file:../kickstart.note::*imports][imports:1]]
use crate::common::*;
use crate::core::*;

use gut::cli::*;
use linefeed::{Interface, ReadResult};
// imports:1 ends here

// [[file:../kickstart.note::659eb3c7][659eb3c7]]
/// A commander for interactive interpreter
#[derive(Default)]
pub struct Command {
    // --
}

impl Command {
    pub fn new() -> Self {
        Self { ..Default::default() }
    }
}

#[derive(StructOpt, Debug)]
pub enum Action {
    /// Quit REPL shell.
    #[structopt(name = "quit", alias = "q", alias = "exit")]
    Quit {},

    /// Show available commands.
    #[structopt(name = "help", alias = "h", alias = "?")]
    Help {},

    /// Terminate program.
    #[structopt(name = "stop", alias = "term", alias = "terminate")]
    Stop {},

    /// Edit population.
    #[structopt(name = "edit", alias = "e")]
    Edit {},
}
// 659eb3c7 ends here

// [[file:../kickstart.note::*core][core:1]]
impl Command {
    fn apply(&mut self, action: &Action) -> Result<()> {
        match action {
            Action::Edit {} => unimplemented!(),
            _ => {
                eprintln!("not implemented yet.");
            }
        }

        Ok(())
    }
}
// core:1 ends here

// [[file:../kickstart.note::86c69857][86c69857]]
/// Return new seeds with user changes in interactive way.
///
/// If return true, user ask for termination.
pub(crate) fn start_repl(control_flag: std::sync::Arc<JobFlag>) -> Result<bool> {
    let version = env!("CARGO_PKG_VERSION");
    let interface = Interface::new("kickstart REPL")?;
    println!("This is kickstart REPL, version {}.", version);
    println!("Enter \"help\" or \"?\" for available commands.");
    println!("Press Ctrl-D or enter \"quit\" or \"q\" to exit REPL shell.");
    println!("Enter \"stop\" to terminate the program.");
    println!("");
    interface.set_prompt("kickstart> ")?;

    let mut command = Command::new();
    while let ReadResult::Input(line) = interface.read_line()? {
        let line = line.trim();
        if !line.is_empty() {
            interface.add_history(line.to_owned());

            let mut args: Vec<_> = line.split_whitespace().collect();
            args.insert(0, "kickstart>");

            match Action::from_iter_safe(&args) {
                // show subcommands
                Ok(Action::Help {}) => {
                    let mut app = Action::clap();
                    app.print_help();
                    println!("");
                }

                Ok(Action::Quit {}) => {
                    let flag = JobType::Run.flag();
                    control_flag.store(flag, std::sync::atomic::Ordering::SeqCst);
                    break;
                }

                Ok(Action::Stop {}) => {
                    info!("User ask for termination.");
                    let flag = JobType::Stop.flag();
                    control_flag.store(flag, std::sync::atomic::Ordering::SeqCst);
                    return Ok(true);
                }

                // edit population
                Ok(Action::Edit {}) => {
                    let flag = JobType::Edit.flag();
                    control_flag.store(flag, std::sync::atomic::Ordering::SeqCst);
                }

                // apply subcommand
                Ok(x) => {
                    if let Err(e) = command.apply(&x) {
                        eprintln!("{:?}", e);
                    }
                }

                // show subcommand usage
                Err(e) => {
                    println!("{e:?}");
                }
            }
        } else {
            println!("");
        }
    }

    Ok(false)
}

fn write_genomes_to_file(genomes: &[MolGenome]) -> Result<()> {
    dbg!("writing file");
    use gchemol::prelude::*;

    let mols: Vec<_> = genomes.iter().map(|m| m.decode()).collect();
    gchemol::io::write("/tmp/test.mol2", &mols)?;

    Ok(())
}

fn load_genomes_from_file() -> Result<Vec<MolGenome>> {
    dbg!("loading file");
    use gchemol::prelude::*;

    let mols: Vec<_> = gchemol::io::read_all("/tmp/test.mol2")?;
    let seeds: Vec<_> = crate::model::compute(mols)?.iter().map(|mp| mp.encode()).collect();

    Ok(seeds)
}
// 86c69857 ends here
