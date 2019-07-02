// src

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*src][src:1]]
/// Simulated annealing.
pub struct Annealer {
    temperature_high: f64,
    temperature_low: f64,
    cooling_rate: f64,
    temperature: Option<f64>,
}

impl Default for Annealer {
    fn default() -> Self {
        Self {
            temperature_high: 5000.0,
            temperature_low: 2000.0,
            cooling_rate: 0.92,
            temperature: None,
        }
    }
}

impl Annealer {
    pub fn new(th: f64, tl: f64) -> Self {
        assert!(th > tl, "temperature_low is high than temperature_high!");

        Self {
            temperature_high: th,
            temperature_low: tl,
            ..Default::default()
        }
    }

    /// Reset tempeature to initial state.
    pub fn reset(&mut self) {
        self.temperature = Some(self.temperature_high);
    }

    /// Return an iterator over temperature.
    pub fn start(&mut self) -> impl Iterator<Item = f64> + '_ {
        std::iter::from_fn(move || {
            let temp: &mut f64 = self.temperature.get_or_insert(self.temperature_high);
            *temp *= self.cooling_rate;

            if *temp <= self.temperature_low {
                None
            } else {
                Some(*temp)
            }
        })
    }
}

// impl Iterator for Annealer {
//     type Item = f64;

//     fn next(&mut self) -> Option<Self::Item> {
//         let temp: &mut f64 = self.temperature.get_or_insert(self.temperature_high);

//         *temp *= self.cooling_rate;
//         // reset temperature if it is too low
//         if *temp <= self.temperature_low {
//             *temp = self.temperature_low
//         }

//         Some(*temp)
//     }
// }

#[test]
fn test_annealer() {
    let mut ann = Annealer::new(500.0, 100.0);

    for x in ann.start() {
        dbg!(x);
    }
}
// src:1 ends here
