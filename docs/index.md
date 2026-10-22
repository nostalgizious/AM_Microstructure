# Microstructure simulation

AMMicrostructure models grain evolution during PBF-LB/M using a three-dimensional cellular automaton. Finite element temperature histories drive nucleation, remelting and competitive growth. Each grain carries a crystallographic orientation, and decentered octahedra determine its advance through the surrounding cells.

The implementation was developed by Ferdinand O. Schnell for the TUM bachelor’s thesis *Simulating Grain Growth in PBF-LB/M with a Cellular Automaton Model to Determine Fracture Criticality* (2024).

## Explore the model

| Topic | Guide |
| --- | --- |
| Thesis contribution, methods and results | [Research](research.md) |
| Cell states, nucleation and octahedral capture | [Cellular automaton model](model-scope.md) |
| Run the two-pass simulation | [Quickstart](quickstart.md) |
| Supply finite element temperature fields | [FEM integration](fem-integration.md) |
| Work with the simulation state | [Simulation API](architecture.md) |
| Inspect grain evolution and thermal histories | [Simulation figures](figures.md) |

## From temperature history to grain structure

![FEM temperatures drive the cellular automaton and produce grain identities and orientations.](assets/thermal-handoff.svg)

Spatial and temporal interpolation transfers FEM temperatures to the cell centres. The CA updates the state of each cell, creates nuclei according to undercooling and advances existing grains into the solidifying region. Grain identity and orientation persist across timesteps, allowing successive thermal cycles to act on the evolving microstructure.

Grain postprocessing provides size, aspect ratio and orientation descriptors for a microstructure-dependent fracture model. The [research guide](research.md) presents the thesis studies, and the [complete thesis](assets/publications/Ferdinand_O_Schnell_BA.pdf) describes the methods and results.
