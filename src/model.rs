// [[file:../kickstart.note::*imports][imports:1]]
use gosh::model::*;
use gchemol::Molecule;

use crate::common::*;
// imports:1 ends here

// [[file:../kickstart.note::4d04f422][4d04f422]]
struct Calculator {
    bunch_mode: bool,
    bbm: BlackBox,
}

impl Calculator {
    /// Construct a calculator using BlackBox model as configured in `dir`.
    fn new(dir: &str) -> Self {
        let bbm = BlackBox::from_dir(dir).expect("bbm failure");
        Self { bunch_mode: false, bbm }
    }

    fn with_bunch_mode(mut self, enabled: bool) -> Self {
        self.bunch_mode = enabled;
        self
    }

    /// Return the calculated results using BlackBox Model.
    fn calculate(&mut self, mols: &[Molecule]) -> Result<Vec<ModelProperties>> {
        let results: Vec<_> = if self.bunch_mode {
            debug!("Calculate {} molecules in bunch mode.", mols.len());
            self.bbm.compute_bunch(mols).context("opt failed in bundle mode")?
        } else {
            // not possible to parallel
            mols.iter()
                .map(|mol| self.bbm.compute(mol).expect("opt mol failure"))
                .collect()
        };

        Ok(results)
    }
}
// 4d04f422 ends here

// [[file:../kickstart.note::b238b72d][b238b72d]]
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
        let calculators: Vec<_> = (0..n).map(|_| Calculator::new(bbm_dir).with_bunch_mode(bunch_mode)).collect();
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

        let nbunch = if self.bunch_mode { self.n_calculators } else { 1 };
        info!("compute {n_mols} molecules with {nbunch} calculators");
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
// b238b72d ends here

// [[file:../kickstart.note::35efc986][35efc986]]
use crate::core::*;

/// setup a runner based on global config and compute a list of molecules.
pub fn compute(mols: Vec<Molecule>) -> Result<Vec<ModelProperties>> {
    debug!("Computing {} molecules ...", mols.len());
    let config = &crate::config::CONFIG;

    let n = config.number_of_calculators;
    let enable_bunch_mode = config.run_in_bunch_mode;
    let mut bbm = BlackBoxModel::from_dir(&config.bbm_dir)?;
    let mut runner = Runner::new(n, &config.bbm_dir, enable_bunch_mode);
    runner.compute(mols)

    // // FIXME: adhoc hacking
    // mols.into_iter()
    //     .map(|mut mol| {
    //         // use builtin optimizer
    //         // let steps = gosh::optim::optimize_geometry_iter(&mut mol, &mut bbm);
    //         // for progress in steps.take(100) {
    //         //     if progress.fmax < 0.1 {
    //         //         break;
    //         //     }
    //         // }
    //         bbm.compute(&mut mol).map(|mut mp| {
    //             mp.set_molecule(mol);
    //             mp
    //         })
    //     })
    //     .collect()
}
// 35efc986 ends here

// [[file:../kickstart.note::*test][test:1]]
#[test]
#[ignore]
fn test_concurrent_calculation() -> Result<()> {
    let bbm_dir = "/share/apps/gulp/opt_reaxff";
    let mut runner = Runner::new(2, bbm_dir, false);

    // jobs
    let molfile = "/home/ybyygu/Incoming/research/kickstart/files/test2/g000.xyz";
    let mols = gchemol::io::read_all(molfile)?;

    let _ = runner.compute(mols)?;

    Ok(())
}
// test:1 ends here
