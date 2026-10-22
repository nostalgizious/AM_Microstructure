# Simulation interface

The cellular automaton advances a persistent grain state using successive cell-temperature fields. A field has shape `(nx, ny, nz)` and contains finite, nonnegative temperatures in kelvin. Keep the same `SimulationState` between steps so the model can use the previous temperature field when evaluating cooling and melting.

```julia
using AMMicrostructure

config = load_demo_config("examples/directional_solidification/config.toml")
p = parameters(config)
T0 = fill(1550.0, p.nx, p.ny, p.nz) # K
state = SimulationState(p, T0; seed=2026)

Tnext = fill(1549.0, p.nx, p.ny, p.nz) # K
counts = step!(state, Tnext; dt_ms=0.02)

state.ids       # Grain ID in each cell; zero means unassigned
state.states    # Cellular-automaton state labels
state.grains    # Grain orientations
state.octa      # Decentered octahedral growth envelopes
state.time_ms   # Elapsed simulation time
```

Coordinates and growth distances use millimetres, time uses milliseconds, and orientation angles use degrees. Julia arrays are column-major, so x varies fastest. The default initialization seeds grains in the bottom layer. The `solid-substrate` initialization extends the seeded grains throughout a field below solidus.

## Grain growth and capture

Each grain carries a crystallographic orientation. Decentered octahedral envelopes attached to solid cells represent its growth front. The envelopes advance according to the local thermal conditions where mushy neighbours are present. An unassigned mushy cell is captured when its centre falls within a candidate grain's rotated envelope and receives that grain's ID and orientation. Captures are applied after candidate traversal.

The cellular automaton uses four state labels: 0 inactive, 1 active solid, 2 mushy or unassigned, and 3 liquid. State and grain assignment are distinct: an assigned grain ID marks microstructure identity, while the state records the thermal category used by the automaton. Assigned fraction is the proportion of cells with a nonzero grain ID. Display colours encode grain IDs.

## Thermal-field handoff

The FEM adapter accepts nodal temperature histories and maps them to CA cell centres. The structured-mesh implementation uses trilinear interpolation in space and linear interpolation between saved times. Its interface converts mesh coordinates from metres to millimetres and time from seconds to milliseconds. See [FEM integration](fem-integration.md) for the input contract and a replay command.
