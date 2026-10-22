# Run a simulation

Install [Julia](https://julialang.org/downloads/) 1.10 or newer, then run commands from the repository root.

## Two-pass grain evolution

```sh
julia --project=. --threads=8,0 examples/fem_replay/run.jl
```

The calculation follows a seeded substrate through two heating and cooling cycles. The second pass remelts part of the structure formed during the first. Saved FEM temperatures are included, so a CalculiX installation is not required.

The default grid contains 384 × 72 × 240 cells with 0.167 µm spacing. To run a smaller case, select another resolution:

```sh
julia --project=. --threads=8,0 examples/fem_replay/run.jl --resolution fine
julia --project=. --threads=8,0 examples/fem_replay/run.jl --resolution quick
```

Results are written to `output/fem_replay/`:

| Output | Contents |
| --- | --- |
| `preview.svg` | Grain sections and their evolution |
| VTK snapshots | Grain IDs, cell states and temperature fields at five stages |
| `history.csv` | Grain assignment and capture through time |
| `slice_history.csv` | Grain and temperature history of a fixed section |

Open the preview in a browser or inspect the VTK fields in a scientific viewer. See [Simulation figures](figures.md) for the two-pass results.

## Directional solidification

```sh
julia --project=. --threads=8,0 examples/directional_solidification/run.jl
```

A seeded bottom layer grows into an initially liquid region under a directional cooling field. The calculation advances a 32 × 6 × 24 cell grid through 100 steps over 2 ms. Results are written to `output/demo/`.

## Supply a temperature field

Use the [simulation API](architecture.md) to pass successive temperature fields to the cellular automaton. The [FEM guide](fem-integration.md) describes nodal temperature transfer, interpolation and units.
