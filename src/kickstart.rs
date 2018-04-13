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

use gchemol::write_as_xyz;
#[test]
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

    write_as_xyz(&all, "/tmp/test.xyz");
}
// ebe7fb0e-fccc-4f07-b951-1dba89782c22 ends here
