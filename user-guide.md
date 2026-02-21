# kickstart user guide

This guide documents the public Linux binaries distributed in GitHub Releases.

## Commands overview

- `kickstart`: run and inspect the structure-search workflow.
- `kickgen`: generate random structures from an input molecule.
- `kickmut`: generate randomly mutated structures.

## `kickgen`

Generate random structures from one input structure.

Usage:

```bash
kickgen [OPTIONS] <INPFILE> <OUTFILE>
```

Arguments:

- `<INPFILE>`: input molecular structure file.
- `<OUTFILE>`: output file for generated structures.

Options:

- `-n <NUMBER_OF_MOLECULES>`: number of generated structures (default: `1`).

Example:

```bash
kickgen input.mol2 output_gen.mol2 -n 10
```

## `kickmut`

Generate random bond-mutation structures from one input structure.

Usage:

```bash
kickmut [OPTIONS] <INPFILE> <OUTFILE>
```

Arguments:

- `<INPFILE>`: input molecular structure file.
- `<OUTFILE>`: output file for mutated structures.

Options:

- `-n <NUMBER_OF_MOLECULES>`: number of generated structures (default: `1`).
- `-m <MUTATION_DEGREE>`: mutation degree per generated structure (default: `1`).

Example:

```bash
kickmut input.mol2 output_mut.mol2 -n 5 -m 2
```

## `kickstart`

Command-line entry for structure exploration and search orchestration.

Usage:

```bash
kickstart [OPTIONS]
```

Options:

- `-p, --print`: print default configuration (TOML).
- `-l, --list`: list calculated items in database.
- `--sort`: sort list output by energy.
- `-r, --run`: run the genetic search workflow.
- `-v, --verbose...`: increase log verbosity.

Examples:

```bash
kickstart --print
kickstart --list
kickstart --list --sort
kickstart --run
```

## Runtime notes

- Official Linux release binaries are static `musl` builds.
- No dynamic runtime dependencies are required on standard x86_64 Linux systems.
