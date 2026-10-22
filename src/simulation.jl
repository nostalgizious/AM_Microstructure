"""
    SimulationState(params, initial_temperature; seed=2026, max_cells=250_000,
                    initialization="bottom-layer")

Initialize the thesis CA with random substrate grains in the bottom cell layer.
`initialization="solid-substrate"` extends those seeded grains through the full
height for a remelting experiment; all initial temperatures must be below solidus.
Temperatures are kelvin, coordinates millimetres, and times milliseconds. A private
RNG and cell-specific capture streams keep seeded results independent of worker count.
"""
mutable struct SimulationState
    params::Params
    ids::Array{Int,3}
    states::Array{Int,3}
    octa::Array{Float64,4}
    grains::Matrix{Float64}
    previous_temperature::Array{Float64,3}
    num_grains::Int
    time_ms::Float64
    step::Int
    rng::Random.Xoshiro
end

function check_temperature(temp, p)
    size(temp) == (p.nx,p.ny,p.nz) || throw(DimensionMismatch("Expected nx × ny × nz temperatures"))
    all(isfinite,temp) || throw(ArgumentError("Temperatures must be finite kelvin values"))
    all(>=(0),temp) || throw(ArgumentError("Kelvin temperatures cannot be negative"))
end

function check_grid(nx,ny,nz,num_grains;max_cells=250_000)
    all(n -> n isa Integer && n > 0,(nx,ny,nz,num_grains,max_cells)) ||
        throw(ArgumentError("Grid sizes, grain count and max_cells must be positive integers"))
    cells = big(nx)*ny*nz
    cells <= max_cells || throw(ArgumentError("Requested $cells CA cells exceeds max_cells=$max_cells. Use a smaller grid; increase max_cells explicitly only after budgeting memory (at least 80 bytes/cell for core arrays, plus working arrays and Julia)."))
    num_grains <= big(nx)*ny || throw(ArgumentError("Too many substrate grains"))
    nothing
end

function SimulationState(p::Params, temp::AbstractArray; seed=2026,max_cells=250_000,initialization="bottom-layer")
    check_grid(p.nx,p.ny,p.nz,p.numGrains;max_cells)
    check_temperature(temp,p)
    initialization in ("bottom-layer","solid-substrate") || throw(ArgumentError("Unknown CA initialization"))
    initialization == "solid-substrate" && maximum(temp) >= p.Tsolidus &&
        throw(ArgumentError("Solid substrate must initially be below solidus"))
    Random.seed!(seed)
    ids,states,octa,grains = initialize(p)
    if initialization == "solid-substrate"
        for k in 2:p.nz
            ids[:,:,k] .= ids[:,:,1]
        end
        states .= 1
    end
    SimulationState(p,ids,states,octa,grains,Float64.(temp),p.numGrains,0.0,0,Random.Xoshiro(seed))
end

"""
    step!(state, temperature; dt_ms=state.params.dt, skip_cold=true)

Advance the original CA with a complete new cell-centred thermal field. Input is
copied, so callers may reuse their buffer. Returns kernel event counts. This is
one-way coupling: latent heat and grain state are not fed back to the FEM solver.
"""
function step!(s::SimulationState, temp::AbstractArray; dt_ms=s.params.dt,skip_cold=true)
    isfinite(dt_ms) && dt_ms > 0 || throw(ArgumentError("dt_ms must be finite and positive"))
    check_temperature(temp,s.params)
    p = s.params
    logvars = Dict{String,Any}()
    if skip_cold && maximum(temp) < p.Tsolidus && maximum(s.previous_temperature) < p.Tsolidus &&
        all(>(0),s.ids) && all(==(1),s.states)
        # No mushy neighbours exist, so the interface-only growth rule freezes
        # all octahedra. Omit microsteps while retaining time, thermal state and RNG.
        logvars["numOctahedralGrowth"] = 0
        logvars["numCellsNucleated"] = 0
        logvars["numInterfaceGrowth"] = 0
        logvars["substeps"] = 0
        logvars["cold_shortcut"] = true
        logvars["numGrains"] = s.num_grains
        logvars["numActiveSolids"] = length(s.ids)
        logvars["currZmin"],logvars["currZmax"],logvars["currHeight"] = 1,p.nz,p.nz
        for key in ("mmmFace","mmpFace","mpmFace","mppFace","pmmFace","pmpFace","ppmFace","pppFace",
                    "numLiquidCells","numCoolingTooQuickly","maxNucleation","timeUpdateCellStates",
                    "timeNucleation","timeUpdateInterfaceCells","timeUpdateGrowingCells","timeToLogVars")
            logvars[key] = 0
        end
    else
        s.num_grains = updateCrystalMeshReduced!(s.ids,s.states,s.octa,s.grains,temp,
            s.previous_temperature,logvars,p,p.nx,p.ny,p.nz,1,p.nz,p.nz,s.num_grains,
            Float64(dt_ms),s.step+1;rng=s.rng)
        logvars["cold_shortcut"] = false
    end
    s.previous_temperature .= temp
    s.time_ms += dt_ms
    s.step += 1
    return logvars
end
