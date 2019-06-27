// config.rs
// :PROPERTIES:
// :header-args: :tangle src/config.rs
// :END:

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*config.rs][config.rs:1]]
use serde::*;
use toml;

#[derive(Debug, Serialize, Deserialize, PartialEq, Clone)]
pub struct Config {
    pub runfile_sp: String,
    pub runfile_opt: String,
    pub molfile: String,
    pub search: Search,
}

#[derive(Debug, Serialize, Deserialize, PartialEq, Clone)]
pub struct Search {
    pub max_generations: u64,
    pub population_size: usize,
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
            },
        }
    }
}

#[test]
fn test_config() {
    let config = Config::default();

    let x = toml::to_string(&config).unwrap();
    println!("{:#?}", x);
}
// config.rs:1 ends here
