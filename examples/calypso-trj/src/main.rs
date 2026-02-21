// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use gchemol::Molecule;
use gut::prelude::*;
use structopt::*;

pub use gut::prelude::*;
// imports:1 ends here

// parse

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*parse][parse:1]]
use gchemol::io;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

struct Data {
    name_of_atom: String,
    number_of_atom: usize,
}

/// read element symbols from input.dat in the parent directory of pso_opt_* file
fn get_symbols_from_pso_input(pso_result_file: &Path) -> Result<HashMap<String, Vec<String>>> {
    let parent = pso_result_file.parent().ok_or(format_err!("wrong path"))?;
    let input_dat = parent.join("../input.dat");

    info!("read CALYPSO input from: {:?}", input_dat.display());
    let txt = read_file(input_dat)?;
    let lines: Vec<_> = txt
        .lines()
        .filter(|line| {
            line.starts_with("NameOfAtoms")
                || line.starts_with("NumberOfAtoms")
                || line.starts_with("NumberOfFormula")
        })
        .take(2)
        .collect();

    assert_eq!(2, lines.len());
    let mut d = HashMap::new();
    for line in lines {
        let parts: Vec<_> = line.split_whitespace().collect();
        let k = &parts[0];
        let v: Vec<_> = parts[2..].iter().map(|v| v.to_string()).collect();
        d.insert(k.to_string(), v);
    }

    Ok(d)
}

/// format calypso result into standard POSCAR file
fn format_pso_opt_file<P: AsRef<Path>>(path: P) -> Result<String> {
    let path = path.as_ref();

    // test if the result file format
    // pso_ini_* contains initial structures
    // pso_opt_* contains optimized structures with final energies
    let filename = path
        .file_name()
        .and_then(|v| v.to_str())
        .ok_or(format_err!("not a normal file"))?;
    let is_ini_file = filename.contains("_ini_");

    if is_ini_file {
        info!("read initial structures from {:?}", path.display());
    } else {
        info!("reading optimized structures from {:?}", path.display());
    }

    let d = get_symbols_from_pso_input(path)?;
    let numbers = &d["NumberOfAtoms"];
    let names = &d["NameOfAtoms"];

    // total number of atoms
    let mut natoms = 0;
    // read in data, get number of atoms in molecule
    let txt = read_file(path)?;
    let nskip = if is_ini_file { 5 } else { 6 };
    if let Some(line) = txt.lines().skip(nskip).take(1).next() {
        natoms = line
            .split_whitespace()
            .map(|s| s.parse::<usize>().expect("number of atoms"))
            .sum();
    } else {
        bail!("incomplete file: {}", path.display());
    }

    // total number of line for one molecule
    let ntotal = if is_ini_file { 7 + natoms } else { 8 + natoms };

    let mut all_lines = txt.lines();

    let mut parts = vec![];
    'out: loop {
        let mut lines = vec![];
        for _ in 0..ntotal {
            if let Some(line) = all_lines.next() {
                lines.push(line.trim_left().to_owned());
            } else {
                break 'out;
            }
        }
        if lines.len() != ntotal {
            break 'out;
        }
        // use energy as the title
        let syms = names.join(" ").to_string();
        if is_ini_file {
            lines.insert(5, syms);
        } else {
            lines.remove(1);
            lines.insert(5, syms);
        }

        let stream = lines.join("\n");
        parts.push(stream);
    }

    // append a blank line
    parts.push("\n".to_string());

    Ok(parts.join("\n"))
}
// parse:1 ends here

// main

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*main][main:1]]
// extern crate gchemol;
// #[macro_use]
// extern crate quicli;

#[derive(Debug, StructOpt)]
struct Cli {
    // Add a positional argument that the user has to supply:
    /// CALYPSO generated result file, e.g. pso_opt_10
    inpfile: String,
    /// The output file name, e.g. /tmp/pso_opt.mol2
    outfile: String,
}

fn main() -> Result<()> {
    let args = Cli::from_args();
    info!("reading from: {}", &args.inpfile);
    let ipath = Path::new(&args.inpfile);
    let opath = Path::new(&args.outfile);
    let stream = format_pso_opt_file(ipath)?;
    write_to_file("/tmp/poscar", &stream)?;

    let mols = io::read("/tmp/poscar")?;

    io::write(&args.outfile, &mols)?;
    info!("wrote to: {}", &args.outfile);

    Ok(())
}
// main:1 ends here
