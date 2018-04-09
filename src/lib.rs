// [[file:~/Workspace/Programming/structure-predication/kick/kickstart.note::0bd0f9c3-bcc6-4163-a8e5-43387b9d65fe][0bd0f9c3-bcc6-4163-a8e5-43387b9d65fe]]
extern crate cgmath;
#[macro_use]
extern crate timeit;
#[macro_use]
extern crate approx;
extern crate rand;
extern crate itertools;

mod kick;
mod ckick;

pub use kick::rand_point_on_sphere;
pub use kick::rand_point_within_sphere;

pub type Point3D = [f64; 3];
pub type Points = Vec<Point3D>;

#[inline]
pub fn euclidean_distance(p1: Point3D, p2: Point3D) -> f64 {
    let mut d2 = 0.0;
    for v in 0..3 {
        let dv = p2[v] - p1[v];
        d2 += dv*dv;
    }

    d2.sqrt()
}
// 0bd0f9c3-bcc6-4163-a8e5-43387b9d65fe ends here
