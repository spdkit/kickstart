// [[file:~/Workspace/Programming/structure-predication/kick/kickstart.note::dacffd90-f296-4155-b1d5-d1b10033e853][dacffd90-f296-4155-b1d5-d1b10033e853]]
type Point3D = [f64; 3];
type Points = Vec<Point3D>;

use rand::{self, Rng};
use rand::distributions::{Distribution, Range};

use itertools::Itertools;

/// create a random point within a sphere
/// References
/// ----------
/// https://stackoverflow.com/a/5408344
fn rand_point_within_sphere(radius: f64) -> Point3D{
    let mut rng = rand::thread_rng();
    let range = Range::new(-radius, radius);

    // using the discarding method, which is simple and also fast
    loop {
        let x = rng.sample(range);
        let y = rng.sample(range);
        let z = rng.sample(range);
        let r2 = x*x + y*y + z*z;
        if r2 <= radius*radius {
            return [x, y, z];
        }
    }
}

/// check if any pair of points come too close
fn close_contact(points: &Points) -> bool {
    let cutoff = 0.4;

    for pair in points.iter().combinations(2) {
        let p1 = pair[0];
        let p2 = pair[1];
        let dx = p2[0] - p1[0];
        let dy = p2[1] - p1[1];
        let dz = p2[2] - p1[2];
        let d2 = dx*dx + dy*dy + dz*dz;
        if d2 <= cutoff {
            println!("{:?}", d2);
            return true
        }
    }

    false
}

fn kick_once(points: &Points, radius: f64) -> Points {
    assert!(radius > 0., "radius should be positive: {}", radius);

    let mut new_points = vec![];
    for &p in points.iter() {
        let d = rand_point_within_sphere(2.);
        let new = [d[0] + p[0], d[1] + p[1], d[2] + p[2]];
        new_points.push(new);
    }

    new_points
}

fn kick(points: Points, radius: f64) -> Points {
    loop {
        let new_points = kick_once(&points, radius);
        if ! close_contact(&new_points) {
            return new_points
        }
    }
}

#[test]
fn test_kick() {
    // positions of atoms in C6H6
    let points = [
        [ -1.25121243e+00,   8.97187180e-01,   0.00000000e+00],
        [  1.43947570e-01,   8.97187180e-01,   0.00000000e+00],
        [  8.41485570e-01,   2.10493818e+00,   0.00000000e+00],
        [  1.43831570e-01,   3.31344718e+00,  -1.19900000e-03],
        [ -1.25099343e+00,   3.31336918e+00,  -1.67800000e-03],
        [ -1.94859443e+00,   2.10516318e+00,  -6.82000000e-04],
        [ -1.80097143e+00,  -5.51298200e-02,   4.50000000e-04],
        [  6.93455570e-01,  -5.53258200e-02,   1.31500000e-03],
        [  1.94116557e+00,   2.10501818e+00,   6.34000000e-04],
        [  6.94031570e-01,   4.26559018e+00,  -1.25800000e-03],
        [ -1.80111543e+00,   4.26565018e+00,  -2.63100000e-03],
        [ -3.04819843e+00,   2.10534618e+00,  -8.62000000e-04]
    ];

    let x = kick(points.to_vec(), 0.8);
    println!("{:?}", x);
}
// dacffd90-f296-4155-b1d5-d1b10033e853 ends here
