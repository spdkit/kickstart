# kickstart

`kickstart` is a command-line toolkit for molecular structure exploration.

This repository is published as a **research artifact** for supplementary information (SI).
The practical distribution path is prebuilt Linux binaries.

## What is included

- `kickstart`: run genetic-search style structure exploration.
- `kickgen`: generate random structures from an input molecule.
- `kickmut`: generate random bond-mutation structures.

## Scope and reproducibility

- This project is provided in a **binary-first** mode.
- Some development dependencies in the source tree are private and are not published.
- For SI verification, use the Linux release binaries attached on GitHub Releases.

## Quick start (Linux x86_64)

1. Download binaries from the latest GitHub Release.
2. Make them executable:

```bash
chmod +x kickstart kickgen kickmut
```

3. Check they run:

```bash
./kickstart --help
./kickgen --help
./kickmut --help
```

## Minimal commands

```bash
./kickgen input.mol2 output_gen.mol2 -n 10
./kickmut input.mol2 output_mut.mol2 -n 5 -m 2
./kickstart --print
```

## Citation

If you use this artifact in academic work, please cite the associated paper and repository release tag used in SI.
