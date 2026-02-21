# kickstart

`kickstart` is a command-line toolkit for stochastic molecular structure generation and refinement.

This repository is published as a **source-available, binary-distributed research artifact** for the Supplementary Information (SI) of the manuscript:

> *Closing the gap between theory and experiment in vibrational spectroscopy via physics-informed theory engineering*

## Scientific context

In the SI workflow (Supplementary Note 1), `kickstart` is used to explore CHON molecular space (typically < 12 heavy atoms) by combining:

- random coordinate kicking,
- distance-geometry filtering,
- evolutionary refinement,
- force-field-based geometry optimization.

This generation stage is part of the Theory Engineering pipeline used to build a quasi-experimental IR library.

## What is included

- `kickstart`: orchestration and search workflow entry point.
- `kickgen`: random structure generation from an input molecule.
- `kickmut`: random bond-mutation structure generation.

## Distribution model and reproducibility

- This project is distributed in **binary-first mode** for Linux x86_64.
- Official release binaries are built with `musl`, statically linked, and have no runtime dynamic-library dependency.
- Source code is provided for transparency and method inspection.
- For SI verification, use the release binaries and the release-specific example assets.

## Documentation

- User guide: `user-guide.md`

## Quick start (Linux x86_64)

1. Download assets from the latest GitHub Release.
2. Make binaries executable:

```bash
chmod +x kickstart kickgen kickmut
```

3. Sanity check:

```bash
./kickstart --help
./kickgen --help
./kickmut --help
```

## Minimal usage

```bash
./kickgen input.mol2 output_gen.mol2 -n 10
./kickmut input.mol2 output_mut.mol2 -n 5 -m 2
./kickstart --print
```

## Citation

If you use this artifact, please cite:

1. the associated manuscript, and
2. the exact GitHub release tag used in your SI or benchmark.
