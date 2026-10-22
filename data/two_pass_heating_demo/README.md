# Retained two-pass thermal history

This dataset contains a CalculiX CrunchiX 2.23 solution over a 64 × 12 × 40 µm specimen. The mesh has 256 DC3D8 elements and 459 nodes. Two stationary prescribed-temperature profiles are centred at x = 24 and 40 µm. Each ramps from a 1500 K baseline to 1800 K over 0.2 ms, holds for 0.2 ms, then cools for 1.6 ms. The sequence spans 4 ms with 200 saved nodal snapshots at 0.02 ms intervals.

Replay from the repository root:

```sh
julia --project=. --threads=8,0 examples/fem_replay/run.jl --resolution quick
```

The run uses this saved thermal history and writes the CA output under `output/fem_replay/`. A new thermal solution can be generated with a working CalculiX executable:

```sh
julia --project=. --threads=8,0 integrations/calculix/run.jl \
  --ccx /path/to/ccx --case two-pass-heating --resolution quick
```

| File | Contents |
| --- | --- |
| `thermal.inp` | Four-step solver input in SI units |
| `thermal.dat` | Raw nodal temperature output in kelvin |
| `solver.log` | Solver completion log |
| `mesh.toml` | Mesh coordinates and node-ID grid |
| `initial.toml` | Initial temperatures and CA seed |
| `config.toml` | Base CA configuration |
| `provenance.toml` | Solver identity and checksums for retained files |

The [FEM guide](../../docs/fem-integration.md) describes the handoff and units. Further solver and provenance details are in [`demo/notes.md`](../../demo/notes.md).
