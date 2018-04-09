// [[file:~/Workspace/Programming/structure-predication/kick/kickstart.note::ebe7fb0e-fccc-4f07-b951-1dba89782c22][ebe7fb0e-fccc-4f07-b951-1dba89782c22]]
use cgmath::prelude::*;
use cgmath::{Quaternion, Vector3};
use rand_point_on_sphere;
use Point3D;
use Points;
use euclidean_distance;

fn random_rotate(points: Vec<Point3D>) -> Points{
    let radius = 1.0;
    let p = rand_point_on_sphere(radius);
    let v = Vector3::from(p);
    let s = (v.magnitude2() + radius*radius).sqrt();
    let rot = Quaternion::from_sv(radius/s, v/s);

    let mut rpoints = vec![];
    for &p in points.iter() {
        let v = Vector3::from(p);
        let t: Point3D = (rot*v).into();
        rpoints.push(t);
    }

    rpoints
}

#[test]
fn test_random_rotate() {
    let points = [[-0.02264019, -0.01300713, -0.06295011],
                  [ 1.37326881, -0.01300713, -0.06295011],
                  [-0.44222819, -0.73391213,  0.82834789],
                  [-0.79257355, -1.33584955, -1.69845937],
                  [-0.76587962,  1.29543401, -0.06295011],
                  [-1.46366314,  1.28242565, -0.77914068],
                  [-1.20324889,  1.43034987,  0.82615384]];

    let rpoints = random_rotate(points.to_vec());

    let npoints = points.len();
    assert_eq!(npoints, rpoints.len());
    assert!(points[0][0] != rpoints[0][0]);

    for i in 0..npoints {
        for j in i..npoints {
            let d1 = euclidean_distance(points[i], points[j]);
            let d2 = euclidean_distance(rpoints[i], rpoints[j]);
            assert_relative_eq!(d1, d2, epsilon=1e-4);
        }
    }
}
// ebe7fb0e-fccc-4f07-b951-1dba89782c22 ends here
