// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use gosh::gchemol;

use gchemol::prelude::*;
use gchemol::Molecule;
use spdkit::random::*;

use crate::common::*;
// imports:1 ends here

// core

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*core][core:1]]
pub(crate) fn mutate_molecule(mol: &mut Molecule) -> Result<Molecule> {
    let nbonds = mol.nbonds();
    if nbonds > 0 {
        // remove bonds randomly to break molecule into parts
        let mut edges: Vec<_> = mol.bonds().map(|b| b.index()).collect();
        let mut rng = thread_rng();
        edges.shuffle(&mut rng);
        let nremoved = nbonds.min(5).max(2);
        let n = rng.gen_range(1, nremoved);
        info!("will remove {} bonds ...", n);
        for i in 0..n {
            mol.remove_bond(edges[i]);
        }
    } else {
        warn!("molecule has no bonds!");
    }

    crate::kick(&mol)
}
// core:1 ends here
