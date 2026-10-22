# Saved CalculiX conduction case

This dataset contains a retained CalculiX CrunchiX 2.23 thermal solution for replay by the cellular automaton. It is a small conduction case with a 2 × 2 × 4 DC3D8 mesh (45 nodes) and 576 CA cells. The saved history covers 2 ms.

Replay from the repository root:

```sh
julia --project=. --threads=1 examples/fem_replay/run.jl --dataset data/calculix_demo
```

The replay maps the saved nodal temperatures to CA cell centres and writes VTK, history, and preview output under `output/fem_replay/`.

| File | Contents |
| --- | --- |
| `thermal.inp` | CalculiX input deck in SI units |
| `thermal.dat` | Nodal temperature history in kelvin |
| `solver.log` | Solver output |
| `mesh.toml` | Mesh axes in metres and node-ID grid |
| `initial.toml` | Initial nodal temperatures and CA initialization settings |
| `config.toml` | CA grid, material parameters, timestep, and random seed |
| `provenance.toml` | Solver and input-file provenance |
| `reference_run.toml` | Recorded replay summary |

The adapter performs the metre-to-millimetre and second-to-millisecond conversions described in the [FEM guide](../../docs/fem-integration.md). Additional dataset and runtime notes are in [`demo/notes.md`](../../demo/notes.md).
