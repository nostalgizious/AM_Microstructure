# Saved layer-heating dataset

This directory contains a retained CalculiX CrunchiX 2.23 thermal solution for solver-free replay. The specimen uses 64 DC3D8 elements and 135 nodes. The 2 ms history begins at 1500 K, heats the top boundary to 1800 K, then releases that boundary for cooling. Material and boundary values are supplied in the solver deck.

Replay from the repository root:

```sh
julia --project=. --threads=1 examples/fem_replay/run.jl --dataset data/layer_heating_demo
```

The CA grid has 3,840 cells (32 × 6 × 20) with 0.002 mm spacing and eight seeded initial grains. Replay maps the nodal history to cell centres and writes results under `output/fem_replay/`.

| File | Contents |
| --- | --- |
| `thermal.inp` | Complete CalculiX input deck in SI units |
| `thermal.dat` | Raw nodal temperature output in kelvin |
| `solver.log` | Solver completion log |
| `mesh.toml` | Mesh coordinates in metres and node-ID ordering |
| `initial.toml` | Initial nodal temperatures and seed |
| `config.toml` | CA configuration and random seed |
| `provenance.toml` | Solver identity and SHA256 checksums |

For replay input units and mapping, see the [FEM guide](../../docs/fem-integration.md). Notes on provenance and regeneration are in [`demo/notes.md`](../../demo/notes.md).
