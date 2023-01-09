// [[file:../kickstart.note::058f52f2][058f52f2]]
use super::*;

use vecfx::*;
// 058f52f2 ends here

// [[file:../kickstart.note::758d935b][758d935b]]
use ncollide3d::math::{Point, Vector};
use ncollide3d::shape::Ball;

#[derive(Debug, Clone)]
pub struct Particle {
    ball: Ball<f64>,
    charge: f64,
    position: Point<f64>,
    velocity: Vector<f64>,
    mass: f64,
}

impl Particle {
    fn new(radius: f64) -> Self {
        Self {
            ball: Ball::new(radius),
            charge: 0.0,
            position: Point::new(0.0, 0.0, 0.0),
            velocity: Vector::new(0.0, 0.0, 0.0),
            mass: 1.0,
        }
    }

    pub fn set_charge(&mut self, charge: f64) {
        self.charge = charge;
    }

    pub fn set_position(&mut self, coord: [f64; 3]) {
        self.position = coord.into();
    }

    pub fn set_mass(&mut self, mass: f64) {
        assert!(mass.is_sign_positive(), "invalid mass: {mass:?}");
        self.mass = mass;
    }

    pub fn set_velocity(&mut self, velocity: [f64; 3]) {
        self.velocity = velocity.into();
    }
}
// 758d935b ends here

// [[file:../kickstart.note::3fb2306d][3fb2306d]]
impl Particle {
    /// Compute particle's contribution to electrostatic energy and force
    ///
    /// # Parameters
    ///   * other: other Particle object
    pub fn electrostatic_energy_and_force(&self, other: &Self) -> (f64, Vector<f64>) {
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
            let (eij, fij) = pi.electrostatic_energy_and_force(pj);
            fi += fij;
            e += eij;
        }
        f.push([fi.x, fi.y, fi.z]);
    }

    let forces = f.as_flat().to_owned();
    (e, forces)
}

/// Update velocity and positions in Velocity Verlet Algorithm
fn velocity_verlet_update<'a, U>(velocity: &mut [f64], pot: &mut Dynamics<'a, U>, dt: f64, masses: &[f64]) -> Result<()> {
    // constant
    let m = masses.as_vector_slice();

    // forces and velocities at time `t`
    let f = pot.get_force()?.as_vector_slice();
    let mut v = velocity.as_vector_slice_mut();

    // Update velocities at t + ∆t/2
    v += 0.5 * dt * f.component_div(&m);

    // Update positions at t + ∆t with velocities
    let dr = &v * dt;

    // Update velocities at t + ∆t
    pot.step_toward(dr.as_slice());
    let f_new = pot.get_force()?.as_vector_slice();
    // NOTE: v is a mutable reference to velocity
    v += 0.5 * dt * f_new.component_div(&m);

    Ok(())
}
// 3fb2306d ends here

// [[file:../kickstart.note::19b64702][19b64702]]
mod charged {
    use super::*;

    #[derive(Debug, Copy, Clone)]
    struct Charged {
        position: Vector<f64>,
        charge: f64,
    }

    impl Charged {
        /// Compute particle's contribution to electrostatic energy and force
        ///
        /// # Parameters
        ///   * other: other Particle object
        pub fn electrostatic_energy_and_force(&self, other: &Self) -> (f64, Vector<f64>) {
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

    fn evaluate_electrostatic_energy_and_force(particles: &[Charged]) -> (f64, Vec<f64>) {
        let n = particles.len();
        let mut e = 0.0;
        let mut f = vec![];
        for i in 0..n {
            let pi = &particles[i];
            let mut fi = Vector::default();
            for j in 0..i {
                let pj = &particles[j];
                let (eij, fij) = pi.electrostatic_energy_and_force(pj);
                fi += fij;
                e += eij;
            }
            f.push([fi.x, fi.y, fi.z]);
        }

        let forces = f.as_flat().to_owned();
        (e, forces)
    }
}
// 19b64702 ends here

// [[file:../kickstart.note::e6da4ab2][e6da4ab2]]
mod collision {
    use super::*;

    use ncollide3d::math::{Isometry, Point, Vector};
    use ncollide3d::query;
    use ncollide3d::shape::{Ball, Compound, Cuboid, ShapeHandle};

    impl Particle {
        /// Compute the time of impact for particle collision within max
        /// allowd time `max_time`. Return None if no collision.
        ///
        /// # Parameters
        ///   * particle: a `Particle` struct
        ///   * max_time: the maximum allowed time of impact.
        pub(super) fn predict_time_of_impact(&self, particle: &Particle, max_time: f64) -> Option<f64> {
            let toi = query::time_of_impact_ball_ball(
                &self.position,     // initial position
                &self.velocity,     // initial velocity
                &self.ball,         // the ball
                &particle.position, // other position
                &particle.velocity, // other velocity
                &particle.ball,     // other ball
                max_time,           // the maximum allowed time of impact.
                0.0,                // target distance between spheres
            )?;

            toi.toi.into()
        }

        /// Update particle velocities after collision for reflection with each other.
        pub(super) fn collision_reflect(&self, other: &Particle) -> (Vector<f64>, Vector<f64>) {
            let ma = self.mass;
            let mb = other.mass;
            let mass_ratio = 1.0 / (1.0 / ma + 1.0 / mb);

            let n = (self.position - other.position).normalize();
            let vi = self.velocity;
            let vo_a = vi - 2.0 * (vi.dot(&n) * n) / ma / mass_ratio;
            let vi = other.velocity;
            let vo_b = vi + 2.0 * (vi.dot(&n) * n) / mb / mass_ratio;

            (vo_a, vo_b)
        }
    }

    fn handle_collision_for_particles(particles: &mut [Particle]) {
        todo!()
    }

    #[test]
    #[ignore]
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
}
// e6da4ab2 ends here

// [[file:../kickstart.note::7fa10d91][7fa10d91]]
use gosh::optim::Dynamics;
use gosh::optim::PotentialOutput;

pub struct System {
    mass: Vec<f64>,
    velocity: Vec<f64>,
}

impl System {
    fn propagate<U>(&mut self, dynamics: &mut Dynamics<U>, dt: f64) -> Result<()> {
        assert_eq!(self.mass.len(), dynamics.position().len());
        velocity_verlet_update(&mut self.velocity, dynamics, dt, &self.mass)?;
        Ok(())
    }

    /// Start MD simulation of particles using with time step `dt` and
    /// max number of iterations `nmax
    pub fn simulate(&mut self, particles: &mut [Particle], dt: f64, nmax: usize) -> Result<()> {
        let mut particles_ = particles.to_vec();
        self.mass = particles.iter().flat_map(|p| [p.mass; 3]).collect();
        self.velocity = particles.iter().flat_map(|p| -> [f64; 3] { p.velocity.into() }).collect();

        let xinit = particles.iter().flat_map(|p| p.position.iter().copied()).collect_vec();
        let mut dynamics = Dynamics::new(xinit.as_slice(), |x: &[f64], f: &mut [f64]| {
            for (i, &p) in x.as_3d().iter().enumerate() {
                particles[i].set_position(p);
            }
            let (energy, forces) = evaluate_energy_and_force(&particles);
            f.clone_from_slice(&forces);
            Ok(energy)
        });

        let np = particles_.len();
        for _ in 0..nmax {
            self.propagate(&mut dynamics, dt)?;
            let positions = dynamics.position().as_3d().iter().copied().collect_vec();
            let velocities = self.velocity.as_3d().iter().copied().collect_vec();
            for i in 0..np {
                particles_[i].set_position(positions[i]);
                particles_[i].set_velocity(velocities[i]);
            }
            for p in (0..np).combinations(2) {
                let (i, j) = (p[0], p[1]);
                let pi = &particles_[i];
                let pj = &particles_[j];
                if let Some(toi) = pi.predict_time_of_impact(pj, dt) {
                    dbg!(toi);
                    let (va, vb) = pi.collision_reflect(pj);
                    particles_[i].set_velocity(va.into());
                    particles_[j].set_velocity(vb.into());
                } else {
                    dbg!("nop");
                }
            }
        }
        // copy velocity back
        drop(dynamics);
        for (v, p) in self.velocity.as_3d().iter().zip(particles) {
            p.set_velocity(*v);
        }

        Ok(())
    }
}
// 7fa10d91 ends here
