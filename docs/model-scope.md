# Cellular automaton model

The model represents microstructure on a three-dimensional cellular-automaton grid. Each cell carries a thermal state and, when assigned, a grain ID. Grain IDs identify discrete grains; colours used in visualizations are categorical labels.

## Cell states

The state labels used by the automaton are:

| Label | Meaning |
| ---: | --- |
| 0 | Inactive |
| 1 | Active solid |
| 2 | Mushy or unassigned |
| 3 | Liquid |

Grain assignment is separate from this state. The assigned fraction is the number of cells with a nonzero grain ID divided by the total cell count.

## Nucleation and growth

Nucleation uses a Gaussian distribution over the change in thermal undercooling. New nuclei receive grain IDs and orientations. Existing grains grow through decentered, rotated octahedral envelopes. A candidate cell is captured when its centre lies inside an eligible envelope; captured cells inherit the grain identity and orientation.

Cell temperatures are supplied at successive timesteps. During cooling, the model evaluates nucleation in the interval between successive temperatures. Growth is active in the solidification region when a solid cell has a mushy-state neighbour. Melting can remove grain assignment, after which cooling can produce new nuclei or regrowth from remaining grains.

## Inputs and interpretation

Temperatures use kelvin, coordinates and growth distances use millimetres, elapsed time uses milliseconds, and orientation angles use degrees. See [Simulation interface](architecture.md) for the API and [FEM thermal handoff](fem-integration.md) for mapped thermal input.

The [simulation figures](figures.md) show the evolution of grain boundaries through successive heating and cooling cycles.
