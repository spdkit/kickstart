// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use gosh::models::*;
use gchemol::Molecule;

use crate::common::*;
// imports:1 ends here

// core

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*core][core:1]]
pub(crate) trait KickModel {
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
