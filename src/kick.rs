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
    let points = [
        [  7.246233,   6.755867,  11.236867],
        [  8.710233,   6.755867,  11.236867],
        [  9.402233,   7.954867,  11.236867],
        [  8.670233,   9.222867,  11.236867],
        [  7.286233,   9.222867,  11.236867],
        [  6.554233,   7.954867,  11.236867],
        [  6.794233,   5.649867,  10.391867],
        [  7.978233,   4.965867,   9.869867],
        [  9.162233,   5.649867,  10.391867],
        [ 10.283233,   5.802867,   9.592867],
        [ 10.587233,   8.116867,  10.391867],
        [  9.402233,  10.167867,  10.391867],
        [  8.710233,  11.061867,   9.592867],
        [  7.246233,  11.061867,   9.592867],
        [  6.554233,  10.167867,  10.391867],
        [  5.369233,   9.483867,   9.869867],
        [  5.369233,   8.116867,  10.391867],
        [  4.941233,   7.069867,   9.592867],
        [  5.673233,   5.801867,   9.592867],
        [  7.978233,   4.471867,   8.575867],
        [  6.794233,   4.632867,   7.730867],
        [  5.673233,   5.279867,   8.224867],
        [  4.941233,   6.224867,   7.379867],
        [  4.489233,   7.330867,   8.224867],
        [  4.489233,   8.624867,   7.730867],
        [  4.941233,   9.730867,   8.575867],
        [  5.673233,  10.675867,   7.730867],
        [  6.793233,  11.322867,   8.224867],
        [ 10.586233,   9.483867,   9.869867],
        [  9.162233,  10.306867,   5.563867],
        [ 10.282233,  10.153867,   6.363867],
        [ 11.014233,   8.885867,   6.363867],
        [ 10.586233,   7.839867,   5.563867],
        [  9.402233,   8.000867,   4.718867],
        [  7.246233,   9.199867,   4.718867],
        [  6.793233,  10.306867,   5.563867],
        [  7.978233,  10.989867,   6.086867],
        [  7.978233,  11.483867,   7.379867],
        [  9.162233,  11.322867,   8.224867],
        [ 10.282233,  10.675867,   7.730867],
        [ 11.467233,   8.624867,   7.730867],
        [ 11.467233,   7.331867,   8.224867],
        [ 11.014233,   6.224867,   7.379867],
        [ 10.586233,   6.471867,   6.085867],
        [  9.402233,   5.787867,   5.563867],
        [  8.670233,   6.732867,   4.718867],
        [  7.285233,   6.732867,   4.718867],
        [  6.553233,   8.000867,   4.718867],
        [  5.673233,  10.153867,   6.363867],
        [  4.941233,   8.885867,   6.363867],
        [  5.369233,   7.839867,   5.563867],
        [  5.369233,   6.471867,   6.086867],
        [  6.554233,   5.787867,   5.563867],
        [  7.246233,   4.894867,   6.363867],
        [  8.710233,   4.894867,   6.363867],
        [  9.162233,   4.632867,   7.730867],
        [ 10.282233,   5.279867,   8.224867],
        [ 11.014233,   9.730867,   8.575867],
        [ 11.014233,   7.069867,   9.592867],
        [  8.710233,   9.199867,   4.718867]];

    let x = kick(points.to_vec(), 1.2);
    println!("{:?}", x);
}
// dacffd90-f296-4155-b1d5-d1b10033e853 ends here
