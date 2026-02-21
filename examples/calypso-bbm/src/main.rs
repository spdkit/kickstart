// [[file:../../../kickstart.note::*imports][imports:1]]
use gut::prelude::*;
use gut::cli::*;

use std::path::PathBuf;
// imports:1 ends here

// [[file:../../../kickstart.note::6be65dd1][6be65dd1]]
pub type Result<T> = std::result::Result<T, Error>;

use gosh::gchemol;
use gosh::model::*;

use gchemol::prelude::*;
use gchemol::Molecule;

fn format_as_gaussian_output(mp: &ModelProperties) -> String {
    let energy = mp.energy.expect("no energy");

    // from eV to Hartree
    let energy = energy / 27.211386024367243;

    let mut lines = String::new();
    let line = format!(" SCF Done:  E(UPBEPBE) = {:20.9}      A.U. after 0 cycles\n", energy);
    lines.push_str(&line);
    lines.push_str("                          Input orientation:\n");
    lines.push_str(" ---------------------------------------------------------------------\n");
    lines.push_str(" Center     Atomic      Atomic             Coordinates (Angstroms)    \n");
    lines.push_str(" Number     Number       Type             X           Y           Z   \n");
    lines.push_str(" ---------------------------------------------------------------------\n");

    // read optimized geometry
    let mol = mp.molecule.as_ref().expect("no molecule");

    let mut i = 1;
    for a in mol.atoms() {
        let [x, y, z] = a.position();
        let line = format!(" {:6}{:>12}{:11}{:16.6}{:12.6}{:12.6}\n", i, a.symbol(), 0, x, y, z);
        lines.push_str(&line);
        i += 1;
    }

    lines.push_str(" ---------------------------------------------------------------------\n");
    lines.push_str("Normal termination of Gaussian 09 at Sat Oct 15 13:07:17 CST 2016.\n");

    lines
}

#[test]
fn test_format_gauss() -> Result<()> {
    let tdir = "/share/apps/gulp/opt_reaxff";
    let mut bbm = BlackBox::from_dir(tdir)?;
    let mol = Molecule::from_file("../../tests/files/c6h6.mol2")?;
    let mp = bbm.compute(&mol)?;
    let s = format_as_gaussian_output(&mp);
    println!("{}", s);

    Ok(())
}
// 6be65dd1 ends here

// [[file:../../../kickstart.note::336598a0][336598a0]]
#[derive(Debug, StructOpt)]
#[clap(author, version, about)]
struct Cli {
    /// Input file: g.xyz
    inpfile: PathBuf,

    /// Output file: gsoutput
    outfile: PathBuf,

    /// BlackBox model profile directory
    bbm_dir: PathBuf,
}

fn main() -> Result<()> {
    let args = Cli::from_args();
    let mol = Molecule::from_file(&args.inpfile)?;

    let mut bbm = BlackBox::from_dir(&args.bbm_dir)?;
    let mp = bbm.compute(&mol)?;
    let s = format_as_gaussian_output(&mp);
    write_to_file(&args.outfile, &s)?;

    Ok(())
}
// 336598a0 ends here
