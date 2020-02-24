// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use gchemol::Molecule;
use gosh::models::*;

use crate::common::*;
// imports:1 ends here

// core

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*core][core:1]]
trait KickModel {
    fn energy(&self) -> f64;
    fn get_energy(&self) -> Result<f64>;
    fn get_optimized_molecule(&self) -> Result<Molecule>;
}

impl KickModel for Molecule {
    fn get_energy(&self) -> Result<f64> {
        let config = &crate::config::CONFIG;
        get_energy(&self, &config.bbm_dir)
    }

    fn energy(&self) -> f64 {
        let config = &crate::config::CONFIG;
        self.get_energy().unwrap_or_else(|e| {
            warn!("found error during energy evaluation:\n {:?}", self.title());
            std::f64::MAX
        })
    }

    fn get_optimized_molecule(&self) -> Result<Molecule> {
        let config = &crate::config::CONFIG;
        get_optimized_molecule(&self, &config.bbm_dir)
    }
}

/// Evaluate energy of molecule in current geometry
fn get_energy(mol: &Molecule, runfile: &str) -> Result<f64> {
    debug!("calculate single point energy using bbm model ...");

    let mut bbm = BlackBox::from_dir(runfile);
    let mr = bbm.compute(mol)?;

    if let Some(energy) = mr.energy {
        debug!("sp energy = {:-12.5}", energy);
        return Ok(energy);
    } else {
        bail!("no energy record found in the output!");
    }
}

/// Return optimized geometry
fn get_optimized_molecule(mol: &Molecule, runfile: &str) -> Result<Molecule> {
    debug!("optimize geometry using bbm model ...");

    // avoid bad geometry
    let mut mol = mol.clone();
    let mut bbm = BlackBox::from_dir(runfile);

    match bbm.compute(&mol) {
        Ok(mr) => {
            if let Some(energy) = mr.energy {
                debug!("opt energy = {:-12.5}", energy);
                if let Some(mol) = mr.molecule {
                    return Ok(mol);
                } else {
                    bail!("no molecule record found in bbm results");
                }
            } else {
                bail!("no energy record found in the output!");
            }
        }

        Err(e) => {
            error!("opt failed, use initial geometry instead.");
            error!("{:?}", e);

            Ok(mol.to_owned())
        }
    }
}
// core:1 ends here

// calculator

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*calculator][calculator:1]]
struct Calculator {
    bunch_mode: bool,
    bbm: BlackBox,
}

impl Calculator {
    /// Construct a calculator using BlackBox model as configured in `dir`.
    fn new(dir: &str) -> Self {
        let bbm = BlackBox::from_dir(dir);
        Self {
            bunch_mode: false,
            bbm,
        }
    }

    fn with_bunch_mode(mut self, enabled: bool) -> Self {
        self.bunch_mode = enabled;
        self
    }

    /// Return the calculated results using Black-Box Model.
    fn calculate(&mut self, mols: &[Molecule]) -> Result<Vec<ModelProperties>> {
        let results: Vec<_> = if self.bunch_mode {
            debug!("Calculate {} molecules in bunch mode.", mols.len());
            self.bbm
                .compute_bunch(mols)
                .with_context(|e| format!("opt failed in bundle mode"))?
        } else {
            // not possible to parallel
            mols.iter()
                .map(|mol| self.bbm.compute(mol).expect("opt mol failure"))
                .collect()
        };

        Ok(results)
    }
}
// calculator:1 ends here

// runner

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*runner][runner:1]]
struct Runner {
    n_calculators: usize,
    bunch_mode: bool,
    runner_tx: Option<crossbeam_channel::Sender<Calculator>>,
    runner_rx: Option<crossbeam_channel::Receiver<Calculator>>,
}

impl Runner {
    /// Construct a runner with `n` calculators.
    fn new(n: usize, bbm_dir: &str, bunch_mode: bool) -> Self {
        use crossbeam_channel::unbounded;

        // put calculators into runner queue
        let (sender, receiver) = unbounded();
        let calculators: Vec<_> = (0..n)
            .map(|_| Calculator::new(bbm_dir).with_bunch_mode(bunch_mode))
            .collect();
        for calc in calculators {
            sender.send(calc).unwrap();
        }

        Self {
            runner_tx: Some(sender),
            runner_rx: Some(receiver),
            n_calculators: n,
            bunch_mode,
        }
    }

    /// Run collected jobs using available calculators.
    fn compute(&mut self, mols: Vec<Molecule>) -> Result<Vec<ModelProperties>> {
        let n_mols = mols.len();
        let sender = self.runner_tx.take().unwrap();
        let receiver = self.runner_rx.take().unwrap();

        let nbunch = if self.bunch_mode {
            self.n_calculators
        } else {
            1
        };
        let bunches: Vec<_> = mols.chunks(nbunch).collect();
        let results: Vec<_> = bunches
            .into_par_iter()
            .map_with((sender, receiver), |(tx, rx), bunch| {
                let mut calculator = rx.recv().unwrap();
                let calculated = calculator.calculate(bunch).unwrap();
                tx.send(calculator).unwrap();
                calculated
            })
            .collect();

        let calculated: Vec<_> = results.concat();
        assert_eq!(n_mols, calculated.len());
        Ok(calculated)
    }
}
// runner:1 ends here

// public

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*public][public:1]]
use crate::core::*;

/// setup a runner based on global config and compute a list of molecules.
pub(crate) fn compute(mols: Vec<Molecule>) -> Result<Vec<ModelProperties>> {
    let config = &crate::config::CONFIG;

    let n = config.number_of_calculators;
    let enable_bunch_mode = config.run_in_bunch_mode;
    let mut runner = Runner::new(n, &config.bbm_dir, enable_bunch_mode);
    runner.compute(mols)
}
// public:1 ends here

// test

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*test][test:1]]
#[test]
#[ignore]
fn test_concurrent_calculation() -> Result<()> {
    let bbm_dir = "/share/apps/gulp/opt_reaxff";
    let mut runner = Runner::new(2, bbm_dir, false);

    // jobs
    let molfile = "/home/ybyygu/Incoming/research/kickstart/files/test2/g000.xyz";
    let mols = gchemol::io::read(molfile)?;

    let _ = runner.compute(mols)?;

    Ok(())
}
// test:1 ends here
