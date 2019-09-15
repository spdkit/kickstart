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
#[derive(Debug, Serialize, Deserialize, PartialEq, Clone)]
pub struct Config {
    pub runfile_sp: String,
    pub runfile_opt: String,
    pub molfile: String,
    pub search: Search,
}

#[derive(Debug, Serialize, Deserialize, PartialEq, Clone)]
pub struct Search {
    pub max_generations: usize,
    pub population_size: usize,
    pub boltzmann_temperature: f64,
    pub mutation_rate: f64,
    pub target_energy: Option<f64>,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            runfile_sp: "/share/apps/mopac/sp".into(),
            runfile_opt: "/share/apps/mopac/opt".into(),
            molfile: "test.mol2".into(),
            search: Search {
                population_size: 10,
                max_generations: 10,
                mutation_rate: 0.1,
                boltzmann_temperature: 10000.0,
                target_energy: None,
            },
        }
    }
}

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
