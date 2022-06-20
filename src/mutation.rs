// [[file:../kickstart.note::*imports][imports:1]]
use gosh::gchemol;

use gchemol::prelude::*;
use gchemol::Molecule;
use spdkit::random::*;

use crate::common::*;
// imports:1 ends here

// [[file:../kickstart.note::7d395c93][7d395c93]]
use gosh::model::ChemicalModel;

fn mutate_by_md(mol: &mut Molecule, pot: &mut impl ChemicalModel) -> Result<()> {
    use surface::docs::dynamics::MolecularDynamics;

    let mut simulation = MolecularDynamics::default();
    simulation.time_step(2.0).temperature(5000.0).softening(true);

    for (i, frame) in simulation.simulate(mol, pot).take(5).enumerate() {
        info!("md step {i:2}, energy = {}", frame?.energy);
    }

    Ok(())
}
// 7d395c93 ends here

// [[file:../kickstart.note::4ebcd3f9][4ebcd3f9]]
pub fn random_bond_mutate(mol: &Molecule, degree: usize) -> Result<Molecule> {
    // import Educate trait here
    use educate::prelude::*;

    let mut new_mol = mol.clone();
    let config = &crate::config::CONFIG;
    for _ in 0..degree {
        if let Some(bbm_dir) = &config.md {
            let mut bbm = gosh::model::BlackBoxModel::from_dir(bbm_dir)?;
            info!("Exploit using minima hopping algorithm with bbm: {bbm_dir:?}");
            mutate_by_md(&mut new_mol, &mut bbm)?;
        }
        new_mol = educate::adhoc::random_bond_mutate(&new_mol);
    }

    Ok(new_mol)
}
// 4ebcd3f9 ends here
