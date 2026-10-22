"""Thesis cellular automaton and reproducible one-way thermal coupling examples."""
module AMMicrostructure

using Random, TOML, Printf, SHA
include("core/crystalGrowth.jl")
include("coupling/arrayOps.jl")
include("demo.jl")
include("simulation.jl")
include("calculix.jl")

export Params, SimulationState, step!, run_demo, load_demo_config, parameters
export ThermalSeries, StructuredThermalMesh, read_calculix_dat, map_temperatures
export write_calculix_case, run_calculix_demo, run_thermal_replay

end
