// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use gosh::models::*;
use gchemol::Molecule;

use crate::common::*;
// imports:1 ends here

// calculator

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*calculator][calculator:1]]
pub(crate) struct Calculator {
    bundled: bool,
    bbm: BlackBox,
}

impl Calculator {
    /// Construct a calculator using BlackBox model as configured in `dir`.
    pub fn new(dir: &str) -> Self {
        let bbm = BlackBox::from_dir(dir);
        Self {
            bundled: false,
            bbm,
        }
    }

    /// Return the calculated results using Black-Box Model.
    pub fn calculate(&mut self, mols: &[Molecule]) -> Result<Vec<ModelProperties>> {
        let results: Vec<_> = if self.bundled {
            self.bbm
                .compute_bunch(mols)
                .with_context(|e| format!("opt failed in bundle mode"))?
        } else {
            // not possible to parallel
            mols.iter()
                .map(|mol| self.bbm.compute(mol).expect("opt mol failed"))
                .collect()
        };

        Ok(results)
    }
}
// calculator:1 ends here

// runner

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*runner][runner:1]]
pub(crate) struct Runner {
    runner_tx: Option<crossbeam_channel::Sender<Calculator>>,
    runner_rx: Option<crossbeam_channel::Receiver<Calculator>>,
}

impl Runner {
    /// Construct a runner with `n` calculators.
    pub(crate) fn new(n: usize, bbm_dir: &str) -> Self {
        use crossbeam_channel::unbounded;

        // put calculators into runner queue
        let (sender, receiver) = unbounded();
        let calculators: Vec<_> = (0..n).map(|_| Calculator::new(bbm_dir)).collect();
        for calc in calculators {
            sender.send(calc).unwrap();
        }

        Self {
            runner_tx: Some(sender),
            runner_rx: Some(receiver),
        }
    }

    /// Run collected jobs using available calculators.
    pub(crate) fn compute(&mut self, mols: Vec<Molecule>) -> Result<Vec<ModelProperties>> {
        let sender = self.runner_tx.take().unwrap();
        let receiver = self.runner_rx.take().unwrap();

        let n = mols.len();
        let results: Vec<_> = mols
            .into_par_iter()
            .map_with((sender, receiver), |(tx, rx), mol| {
                let mut calculator = rx.recv().unwrap();
                let calculated = calculator.calculate(&vec![mol]).unwrap();
                tx.send(calculator).unwrap();
                calculated
            })
            .collect();

        let calculated: Vec<_> = results.concat();
        assert_eq!(n, calculated.len());
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

    // FIXME: config optin for number of concurrent runners
    let mut runner = Runner::new(2, &config.bbm_dir);
    runner.compute(mols)
}
// public:1 ends here

// test

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*test][test:1]]
#[test]
fn test_concurrent_calculation() -> Result<()> {
    let bbm_dir = "/share/apps/gulp/opt_reaxff";
    let mut runner = Runner::new(2, bbm_dir);

    // jobs
    let molfile = "/home/ybyygu/Incoming/research/kickstart/files/test2/g000.xyz";
    let mols = gchemol::io::read(molfile)?;

    let _ = runner.compute(mols)?;

    Ok(())
}
// test:1 ends here
