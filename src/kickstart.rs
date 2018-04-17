// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::ebe7fb0e-fccc-4f07-b951-1dba89782c22][ebe7fb0e-fccc-4f07-b951-1dba89782c22]]
use Point3D;
use Points;
use geometry::rand_point_on_sphere;
use geometry::rand_point_within_sphere;
use geometry::euclidean_distance;
use geometry::rand_rotate;

/// a subgroup of points in parent point clound
#[derive(Clone, Debug)]
struct Fragment<'a> {
    points: &'a Points,
    /// indices to point clound
    members: Vec<usize>,

    /// current geometric center
    center: Point3D,
}

impl <'a> Fragment<'a> {
    fn new(points: &'a Points) -> Self {
        Fragment {
            points: points,
            center: [0.0; 3],
            members: vec![],
        }
    }

    /// get positions of points in fragment
    fn positions(&self) -> Points {
        let pts: Points = self.members.iter().map(|&x| self.points[x]).collect();

        pts
    }

    /// kick the fragment to a new location randomly
    /// Parameters
    /// ----------
    /// radius: the sphere radius
    ///
    fn kick(&self, radius: f64) -> Points {
        let pts = rand_rotate(self.positions());

        let center = rand_point_within_sphere(radius);
        let pts = translated(&pts, center);

        pts
    }
}

/// initialize points and members
impl <'a> From<&'a Points> for Fragment<'a> {
    fn from(points: &'a Points) -> Self {
        let indices: Vec<usize> = (0..points.len()).collect();
        let mut frag = Fragment::new(points);
        frag.members = indices;

        frag
    }
}

/// get positions when translated to new location
fn translated(points: &Points, loc: Point3D) -> Points {
    let mut pts = vec![];
    for &pt in points.iter() {
        let mut p = [0.0; 3];
        for v in 0..3 {
            p[v] = pt[v] + loc[v];
        }
        pts.push(p);
    }

    pts
}

#[test]
#[ignore]
fn test_fragment() {
    let points = [[-0.02264019, -0.01300713, -0.06295011],
                  [ 1.37326881, -0.01300713, -0.06295011],
                  [-0.44222819, -0.73391213,  0.82834789],
                  [-0.79257355, -1.33584955, -1.69845937],
                  [-0.76587962,  1.29543401, -0.06295011],
                  [-1.46366314,  1.28242565, -0.77914068],
                  [-1.20324889,  1.43034987,  0.82615384]];

    let points = points.to_vec();

    let max_radius = points.len() as f64 * 2.0 * 4.0;

    let frag = Fragment::from(&points);
    let mut all = vec![];
    for _ in 0..500 {
        let positions = frag.kick(max_radius);
        for &pt in positions.iter() {
            all.push(pt);
        }
    }
}
// ebe7fb0e-fccc-4f07-b951-1dba89782c22 ends here

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::541c7076-3722-47ce-a867-9bb88122e954][541c7076-3722-47ce-a867-9bb88122e954]]
use nalgebra::{self, Isometry3, Vector3};
use ncollide::shape::{Compound3, Ball, Cuboid, ShapeHandle};

use gchemol::{
    Molecule,
};

// convert molecule to a 3D shape compound
fn molecule_to_3dshape(mol: &Molecule) -> Compound3<f64> {
    let delta = |position: Point3D| {
        Isometry3::new(Vector3::new(position[0],
                                    position[1],
                                    position[2]),
                       nalgebra::zero())
    };

    // 1) Initialize the shape list.
    let mut shapes = Vec::new();
    for a in mol.atoms() {
        println!("{:?}", a);
        let d = delta(a.position);
        let ball = ShapeHandle::new(Ball::new(1.0));
        shapes.push((d, ball));
    }

    // 2) Create the compound shape.
    Compound3::new(shapes)
}
// 541c7076-3722-47ce-a867-9bb88122e954 ends here

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::3b01d553-3e90-4a1d-ade5-ab3f71030849][3b01d553-3e90-4a1d-ade5-ab3f71030849]]
use rand::{thread_rng, Rng};
use petgraph::algo::has_path_connecting;
use std::collections::HashMap;

const EPSILON: f64 = 1.0E-6;

fn get_bound_matrix(mol: &Molecule) -> Vec<Vec<f64>>{
    let mut dm = mol.distance_matrix();
    println!("{:?}", dm.len());
    let max_rij = 90.0;
    for i in mol.graph.node_indices() {
        for j in mol.graph.node_indices() {
            if i < j {
                if let Some(e) = mol.graph.find_edge(i, j) {
                    //
                } else {
                    let ai = &mol.graph[i];
                    let aj = &mol.graph[j];
                    let ri = ai.vdw_radius().unwrap();
                    let rj = aj.vdw_radius().unwrap();
                    let rij = ri + rj;
                    if ! has_path_connecting(&mol.graph, i, j, None) {
                        let i = i.index();
                        let j = j.index();
                        dm[i][j] = rij;
                        // if dm[j][i] < rij {
                        //     dm[j][i] = 5.0*rij;
                        // }
                        dm[j][i] = max_rij;
                        // dm[j][i] += rij;
                    }
                }
            }
        }
    }

    dm
}

fn clean_structures(mols: &mut Vec<Molecule>, bounds: &Vec<Vec<f64>>) {
    let mut rng = thread_rng();
    let indices: Vec<_> = (0..mols.len()).collect();

    let mut all = vec![];
    for ref mol in mols.iter() {
        let atoms: Vec<_> = mol.atoms().map(|a| a.index).collect();
        all.push(atoms);
    }

    let mut rng = thread_rng();
    let maxlam = 1.0;
    let minlam = 0.001;
    let maxcycle = 1000;
    let maxstep = 500;
    let dlam = (maxlam - minlam)/(maxcycle as f64);
    let mut lam = maxlam;
    // enter spe loop
    for icycle in 0..maxcycle {
        lam -= dlam;
        // choose two different molecule
        let (i, j) = loop {
            let &i = rng.choose(&indices).unwrap();
            // choose atom index
            let &j = rng.choose(&indices).unwrap();
            if i != j {
                break (i,j);
            }
        };

        // choose atom index from molecule i
        let &ni = rng.choose(&all[i]).expect("ni");
        let &pi = &mols[i].get_atom(ni).expect("pi").position;
        for istep in 0..maxstep {
            println!("cycle {:}, step {:}", icycle, istep);
            // choose atom index from molecule j
            let &nj = rng.choose(&all[j]).expect("nj");
            let &pj = &mols[j].get_atom(nj).expect("pj").position;
            let dij = euclidean_distance(pi, pj);

            let bij = bounds[ni.index()][nj.index()];
            let bji = bounds[nj.index()][ni.index()];
            let mut bound = [bij, bji];
            if bij > bji {
                bound.swap(0, 1);
            }

            // update coordinates
            let lij = bound[0];
            let uij = bound[1];
            if dij >= lij && dij <= uij {
                continue;
            }

            println!("{:?}", (i,j, ni, nj, dij, bound));
            println!("{:?}", (pi, pj));
            // let rij = if (dij - lij).abs() < (dij - uij).abs() {
            //     uij
            // } else {
            //     lij
            // };
            let rij = lij;

            // calculate position shift
            let mut si = [0.0; 3];
            let mut sj = [0.0; 3];
            let disparity = lam*(rij - dij) / (dij + EPSILON);
            for v in 0..3 {
                si[v] += 0.5 * disparity * (pi[v] - pj[v]);
                sj[v] += 0.5 * disparity * (pj[v] - pi[v]);
            }

            // translate all atom in molecule i
            {
                for &vi in all[i].iter() {
                    let mut atom = mols[i].get_atom_mut(vi).expect("i, vi");
                    for v in 0..3 {
                        atom.position[v] += si[v];
                    }
                }
            }
            // translate all atom in molecule j
            {
                for &vj in all[j].iter() {
                    let mut atom = mols[j].get_atom_mut(vj).expect("j, vj");
                    for v in 0..3 {
                        atom.position[v] += sj[v];
                    }
                }
            }
        }
    }
}

#[test]
fn test_clean_structures() {
    let mut mol = Molecule::from_file("/tmp/test.mol2").unwrap();
    let bounds = get_bound_matrix(&mol);
    let mut mols = mol.fragment();
    clean_structures(&mut mols, &bounds);
    let mut atoms = HashMap::new();
    for ref mol in mols.iter() {
        for a in mol.atoms() {
            atoms.insert(a.index, a.clone());
        }
    }
    let mut keys: Vec<_> = atoms.keys().collect();
    keys.sort();
    let mut newmol = Molecule::new("final");
    for n in keys {
        newmol.add_atom(atoms[n].clone());
    }
    newmol.to_file("/tmp/test3.mol2");
    println!("{:?}", bounds);
}
// 3b01d553-3e90-4a1d-ade5-ab3f71030849 ends here

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::6af98c4f-b799-46ae-a70d-09c54b5c1283][6af98c4f-b799-46ae-a70d-09c54b5c1283]]
#[test]
fn test_molecule_to3dshape(){
    use ncollide::query;

    let mol = Molecule::from_file("/home/ybyygu/Workspace/Programming/gchemol/tests/data/c2h4.xyz").unwrap();
    let s1 = molecule_to_3dshape(&mol);
    let s2 = molecule_to_3dshape(&mol);

    // initialize positions
    let p1 = Isometry3::new(Vector3::new(3.0, 3.0, 3.0), nalgebra::zero());
    let p2 = Isometry3::new(Vector3::new(1.0, 1.0, 1.0), nalgebra::zero());

    // initialize velocities
    let v1 = Vector3::new(2.0, 2.0, 2.0);
    let v2 = Vector3::new(-0.5, -0.5, -0.5);

    let toi = query::time_of_impact(
        &p1,
        &v1,
        &s1,
        &p2,
        &v2,
        &s2,
    );

    println!("{:?}", toi);
}
// 6af98c4f-b799-46ae-a70d-09c54b5c1283 ends here

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::07e66245-e4ae-4f4f-9565-a8cccfb88c56][07e66245-e4ae-4f4f-9565-a8cccfb88c56]]
#[test]
#[ignore]
fn test_molecule_spe() {
    use gchemol::Molecule;
    use spe::pspe;

    let mut mol = Molecule::from_file("/tmp/test.mol2").unwrap();
    for b in mol.bonds() {
        println!("{:?}", b);
    }

    let mut dm = mol.distance_matrix();
    println!("{:?}", dm.len());
    let max_rij = 40.0;
    for i in mol.graph.node_indices() {
        for j in mol.graph.node_indices() {
            if i < j {
                if let Some(e) = mol.graph.find_edge(i, j) {
                    //
                } else {
                    let ai = &mol.graph[i];
                    let aj = &mol.graph[j];
                    let ri = ai.vdw_radius().unwrap();
                    let rj = aj.vdw_radius().unwrap();
                    let rij = ri + rj;
                    if ! has_path_connecting(&mol.graph, i, j, None) {
                        let i = i.index();
                        let j = j.index();
                        dm[i][j] = rij;
                        // if dm[j][i] < rij {
                        //     dm[j][i] = 5.0*rij;
                        // }
                        dm[j][i] += rij;
                    }
                }
            }
        }
    }

    let mut positions: Vec<_> = mol.positions().map(|v| *v).collect();

    let maxcycle = 500;
    let maxstep = 1000;
    let maxlam = 1.0;
    let minlam = 0.01;
    pspe(&mut positions, &dm, maxcycle, maxstep, maxlam, minlam);
    mol.set_positions(positions);
    mol.to_file("/tmp/test2.mol2");
    println!("{:?}", dm);
}
// 07e66245-e4ae-4f4f-9565-a8cccfb88c56 ends here
