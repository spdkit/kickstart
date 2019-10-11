// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use serde::*;
use toml;
// imports:1 ends here

// global

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*global][global:1]]
lazy_static! {
    /// Global settings.
    pub static ref CONFIG: Config = {
        let config_file = format!("{}.conf", env!("CARGO_PKG_NAME"));
        println!("configfile {}", config_file);

        let toml_str = quicli::fs::read_file(config_file).expect("Failed to read config file!");
        toml::from_str(&toml_str).expect("Failed to parse toml config!")
    };
}
// global:1 ends here

// base

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*base][base:1]]
#[derive(Debug, Serialize, Deserialize)]
pub struct Config {
    /// The path to BlackBox Model (bbm) directory.
    pub bbm_dir: String,
    /// The number of calculators run in parallel
    pub number_of_calculators: usize,
    /// If run in bunch mode, a bunch of molecules will be put into a single
    /// calculation, useful for cheap method such as force field.
    pub run_in_bunch_mode: bool,
    /// The path to a file containing initial molecule with multiple fragments
    /// (based on connectivity)
    pub molfile: String,
    /// Evolution search parameters
    pub search: Search,
}

#[derive(Debug, Serialize, Deserialize, Default)]
pub struct Search {
    pub max_generations: usize,
    pub population_size: usize,
    pub boltzmann_temperature: f64,
    pub mutation_rate: f64,
    pub target_energy: Option<f64>,
    /// the last geenrations for termination test
    pub termination_nlast: Option<usize>,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            bbm_dir: "/share/apps/mopac/opt".into(),
            molfile: "test.mol2".into(),
            number_of_calculators: 1,
            run_in_bunch_mode: false,
            search: Search {
                population_size: 10,
                max_generations: 10,
                mutation_rate: 0.1,
                boltzmann_temperature: 10000.0,
                ..Default::default()
            },
        }
    }
}

/// Print default configuration.
pub fn print_default_config() {
    let x = toml::to_string(&Config::default()).unwrap();
    println!("{:}", x);
}

#[test]
fn test_config() {
    let config = Config::default();

    let x = toml::to_string(&config).unwrap();
    println!("{:#?}", x);
}
// base:1 ends here
