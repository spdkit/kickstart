// [[file:../kickstart.note::058f52f2][058f52f2]]
use super::*;
// 058f52f2 ends here

// [[file:../kickstart.note::758d935b][758d935b]]
pub struct Particle {
    ball: Ball<f64>,
    charge: f64,
    position: Point<f64>,
    velocity: Vector<f64>,
}

impl Particle {
    fn new(radius: f64) -> Self {
        Self {
            ball: Ball::new(radius),
            charge: 0.0,
            position: Point::new(0.0, 0.0, 0.0),
            velocity: Vector::new(0.0, 0.0, 0.0),
        }
    }

    pub fn set_charge(&mut self, charge: f64) {
        self.charge = charge;
    }

    pub fn set_position(&mut self, coord: [f64; 3]) {
        self.position = coord.into();
    }

    pub fn set_velocity(&mut self, velocity: [f64; 3]) {
        self.velocity = velocity.into();
    }

    /// 计算该点与其它点之间的静电相互作用力.
    ///
    /// # Parameters
    ///   * other: other Particle object
    pub fn elastic_energy_and_force(&self, other: &Self) -> (f64, Vector<f64>) {
        let ke = 8.987551792314E-9;
        let ri = self.position - other.position;
        let q = self.charge;
        let qi = other.charge;
        let ri_norm2 = ri.norm_squared();
        let pre_factor = ke * q * qi;
        let fij = pre_factor / ri_norm2 * ri.normalize();
        let eij = pre_factor / ri_norm2.sqrt();
        (eij, fij)
    }
}
// 758d935b ends here

// [[file:../kickstart.note::7fa10d91][7fa10d91]]
use gosh::optim::Dynamics;
use vecfx::*;

fn evaluate_energy_and_force(particles: &[Particle]) -> (f64, Vec<f64>) {
    use vecfx::*;

    let n = particles.len();
    let mut e = 0.0;
    let mut f = vec![];
    for i in 0..n {
        let pi = &particles[i];
        let mut fi = Vector::default();
        for j in 0..i {
            let pj = &particles[j];
            let (eij, fij) = pi.elastic_energy_and_force(pj);
            fi += fij;
            e += eij;
        }
        f.push([fi.x, fi.y, fi.z]);
    }

    let forces = f.as_flat().to_owned();
    (e, forces)
}

fn simulate() -> Result<()> {
    use gosh::optim::PotentialOutput;
    let p1 = Particle::new(1.0);
    let p2 = Particle::new(1.0);
    let p3 = Particle::new(1.0);
    let mut particles = vec![p1, p2, p3];
    let xinit = particles.iter().flat_map(|p| p.position.iter().copied()).collect_vec();
    let mut dynamics = Dynamics::new(xinit.as_slice(), |x: &[f64], f: &mut [f64]| {
        for (i, &p) in x.as_3d().iter().enumerate() {
            particles[i].set_position(p);
        }
        let (energy, forces) = evaluate_energy_and_force(&particles);
        f.clone_from_slice(&forces);
        Ok(energy)
    });
    todo!()
}
// 7fa10d91 ends here

// [[file:../kickstart.note::3fb2306d][3fb2306d]]
use ncollide3d::math::{Isometry, Point, Vector};
use ncollide3d::query;
use ncollide3d::shape::{Ball, Compound, Cuboid, ShapeHandle};

#[test]
fn test_nc() {
    let ball1 = Ball::new(1.0);
    let ball2 = Ball::new(1.0);

    // 初始化位置
    let ball_pos1 = Point::new(0.0, 0.0, 0.0);
    let ball_pos2 = Point::new(3.0, 3.0, 3.0);

    // 初始化速度
    let ball_vel1 = Vector::new(2.0, 2.0, 2.0);
    let ball_vel2 = Vector::new(-0.5, -0.5, -0.5);

    if let Some(toi) = query::time_of_impact_ball_ball(
        &ball_pos1, // 初始位置
        &ball_vel1, // 初始速度
        &ball1,     // 形状
        &ball_pos2, //
        &ball_vel2, //
        &ball2,     //
        3.0,        // the maximum allowed time of impact.
        0.0,        // target distance between spheres
    ) {
        // 如果是Some(0.0), 表示已经接触
        // 如果是None, 表示不会接触
        // 如果是Some(2.0), 表示2.0单位的时间后会接触.
        println!("{:?}", toi);
    }
}

impl Particle {
    pub fn predict_time_of_impact(&self, particle: &Particle) -> Option<f64> {
        let toi = query::time_of_impact_ball_ball(
            &self.position,     // 初始位置
            &self.velocity,     // 初始速度
            &self.ball,         // 形状
            &particle.position, //
            &particle.velocity, //
            &particle.ball,     //
            3.0,                // the maximum allowed time of impact.
            0.0,                // target distance between spheres
        )?;

        {
            // 如果是Some(0.0), 表示已经接触
            // 如果是None, 表示不会接触
            // 如果是Some(2.0), 表示2.0单位的时间后会接触.
            println!("{:?}", toi);
        }
        toi.toi.into()
    }
}
// 3fb2306d ends here
