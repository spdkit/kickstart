// [[file:~/Workspace/Programming/structure-predication/kick/kickstart.note::dacffd90-f296-4155-b1d5-d1b10033e853][dacffd90-f296-4155-b1d5-d1b10033e853]]
use rand::{self, Rng};
use rand::distributions::{Distribution, Range, Normal};

use Point3D;
use Points;

use itertools::Itertools;

/// create a random point within a sphere
/// References
/// ----------
/// https://stackoverflow.com/a/5408344
pub fn rand_point_within_sphere(radius: f64) -> Point3D {
    let mut rng = rand::thread_rng();
    let range = Range::new(-radius, radius);

    // using the discarding method, which is simple and also fast
    let radius2 = radius*radius;
    loop {
        let x = rng.sample(range);
        let y = rng.sample(range);
        let z = rng.sample(range);
        let r2 = x*x + y*y + z*z;
        if r2 <= radius2 {
            return [x, y, z];
        }
    }
}

/// Generating uniformly distributed point on a sphere
/// Alternative method 1 as described in:
/// http://corysimon.github.io/articles/uniformdistn-on-sphere/
pub fn rand_point_on_sphere(radius: f64) -> Point3D {
    debug_assert!(radius > 0.0, "sphere radius cannot be negative: {:?}", radius);

    let mut rng = rand::thread_rng();

    let radius2 = radius*radius;
    let normal = Normal::new(0.0, 10.0);
    // avoid floating point precision lost when dividing by a very small number.
    let min_radius2 = 0.1;
    loop {
        let x = rng.sample(normal);
        let y = rng.sample(normal);
        let z = rng.sample(normal);

        let r2: f64 = x*x + y*y + z*z;
        if r2 > min_radius2 {
            let s = radius/r2.sqrt();
            return [x*s, y*s, z*s];
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

// #[test]
// fn test_rand_point_on_sphere() {
//     use std::io::prelude::*;
//     use std::fs::File;

//     let mut buffer = File::create("foo.xyz").unwrap();

//     // let mut points = vec![];
//     buffer.write(b"600\n").unwrap();
//     buffer.write(b"title\n").unwrap();
//     for i in 0..600 {
//         let p = rand_point_on_sphere(20.);
//         let s = format!("C {:-12.4}{:-12.4}{:-12.4}\n", p[0], p[1], p[2]);
//         buffer.write(s.as_bytes());
//     }
// }
// dacffd90-f296-4155-b1d5-d1b10033e853 ends here
