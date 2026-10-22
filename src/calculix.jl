# Deliberately bounded adapter: rectangular, axis-aligned linear hexahedra.
# Solver geometry/time use m/s; CA geometry/time use mm/ms; both use K.
struct StructuredThermalMesh
    x_m::Vector{Float64}
    y_m::Vector{Float64}
    z_m::Vector{Float64}
    node_ids::Array{Int,3}
    function StructuredThermalMesh(x,y,z,ids)
        for axis in (x,y,z)
            length(axis) >= 2 && all(isfinite,axis) && all(>(0),diff(axis)) ||
                throw(ArgumentError("Mesh axes must be finite and strictly increasing"))
        end
        size(ids) == (length(x),length(y),length(z)) || throw(DimensionMismatch("Node ID grid"))
        length(unique(ids)) == length(ids) && all(>(0),ids) || throw(ArgumentError("Unique positive node IDs required"))
        new(Float64.(x),Float64.(y),Float64.(z),Int.(ids))
    end
end

"""Nodal temperature snapshots: seconds, explicit node IDs, and kelvin (nodes × times)."""
struct ThermalSeries
    times_s::Vector{Float64}
    node_ids::Vector{Int}
    temperature_K::Matrix{Float64}
    function ThermalSeries(t,ids,temp)
        !isempty(t) && all(isfinite,t) && all(>=(0),t) && all(>(0),diff(t)) ||
            throw(ArgumentError("Snapshot times must be finite, nonnegative and strictly increasing"))
        !isempty(ids) && all(>(0),ids) && length(unique(ids)) == length(ids) ||
            throw(ArgumentError("Unique positive node IDs required"))
        size(temp) == (length(ids),length(t)) || throw(DimensionMismatch("Temperature matrix must be nodes × times"))
        all(isfinite,temp) && all(>=(0),temp) || throw(ArgumentError("Invalid kelvin temperature"))
        new(Float64.(t),Int.(ids),Float64.(temp))
    end
end

parse_ccx_float(s) = parse(Float64,replace(s,'D'=>'E','d'=>'e'))

"""
    read_calculix_dat(path; node_set="NALL")

Read scalar temperature blocks emitted by `*NODE PRINT,NSET=NALL` / `NT`.
Requires exactly one full block for each strictly increasing total time. Other
output blocks and node sets are ignored. Missing/duplicate nodes are errors.
This is not an FRD parser or a general CalculiX mesh importer.
"""
function read_calculix_dat(path; node_set="NALL")
    times = Float64[]
    blocks = Dict{Int,Float64}[]
    active = false
    for line in eachline(path)
        header = match(r"^\s*temperatures(?:\s*\([^)]*\))?\s+for set\s+(\S+)\s+and time\s+(\S+)"i,line)
        if header !== nothing
            active = uppercase(header[1]) == uppercase(node_set)
            if active
                push!(times,parse_ccx_float(header[2]))
                push!(blocks,Dict{Int,Float64}())
            end
            continue
        end
        active || continue
        fields = split(strip(line))
        isempty(fields) && continue
        if length(fields) == 2 && tryparse(Int,fields[1]) !== nothing
            id = parse(Int,fields[1])
            haskey(last(blocks),id) && error("Duplicate node $id in temperature block")
            last(blocks)[id] = parse_ccx_float(fields[2])
        else
            active = false
        end
    end
    isempty(blocks) && error("No temperature output for node set $node_set in $path")
    ids = sort!(collect(keys(first(blocks))))
    isempty(ids) && error("Empty temperature block")
    all(sort!(collect(keys(b))) == ids for b in blocks) || error("Inconsistent/missing temperature nodes")
    ThermalSeries(times,ids,[b[id] for id in ids, b in blocks])
end

function nodal_temperatures(series::ThermalSeries,time_s)
    isfinite(time_s) || throw(ArgumentError("Time must be finite"))
    lo,hi = first(series.times_s),last(series.times_s)
    tol = 32eps(max(abs(lo),abs(hi),1.0))
    lo-tol <= time_s <= hi+tol || throw(ArgumentError("Requested time lies outside available snapshots"))
    t = clamp(time_s,lo,hi)
    left = searchsortedlast(series.times_s,t)
    left == length(series.times_s) && return copy(series.temperature_K[:,left])
    fraction = (t-series.times_s[left])/(series.times_s[left+1]-series.times_s[left])
    (1-fraction).*series.temperature_K[:,left] .+ fraction.*series.temperature_K[:,left+1]
end

function axis_bracket(axis,value)
    first(axis) <= value <= last(axis) || throw(ArgumentError("CA cell lies outside FEM domain"))
    left = min(searchsortedlast(axis,value),length(axis)-1)
    left,(value-axis[left])/(axis[left+1]-axis[left])
end

"""
    map_temperatures(mesh, series, params, time_ms)

Trilinearly interpolate nodal temperatures to CA cell centres `(i-0.5)*dx`.
Linearly interpolate between solver snapshots in time. Reject extrapolation and
missing nodes; the CA origin is (0,0,0). Only axis-aligned rectangular grids are
supported. For general tetrahedra use barycentric mapping in a separate adapter.
"""
function map_temperatures(mesh::StructuredThermalMesh,series::ThermalSeries,p::Params,time_ms)
    vals = nodal_temperatures(series,time_ms/1000)
    positions = Dict(id=>i for (i,id) in enumerate(series.node_ids))
    all(haskey(positions,id) for id in mesh.node_ids) || throw(ArgumentError("Temperature series misses mesh nodes"))
    nodal = reshape([vals[positions[id]] for id in mesh.node_ids],size(mesh.node_ids))
    x_brackets = [axis_bracket(mesh.x_m,(i-0.5)*p.dx/1000) for i in 1:p.nx]
    y_brackets = [axis_bracket(mesh.y_m,(j-0.5)*p.dx/1000) for j in 1:p.ny]
    z_brackets = [axis_bracket(mesh.z_m,(k-0.5)*p.dx/1000) for k in 1:p.nz]
    mapped = zeros(p.nx,p.ny,p.nz)
    Threads.@threads :static for k in 1:p.nz
      for j in 1:p.ny,i in 1:p.nx
        a,u = x_brackets[i]
        b,v = y_brackets[j]
        c,w = z_brackets[k]
        for dk in 0:1,dj in 0:1,di in 0:1
            weight = (di==0 ? 1-u : u)*(dj==0 ? 1-v : v)*(dk==0 ? 1-w : w)
            mapped[i,j,k] += weight*nodal[a+di,b+dj,c+dk]
        end
      end
    end
    mapped
end

"""Retain only two spatially mapped FEM snapshots and interpolate them in time."""
mutable struct ThermalFieldCache
    mesh::StructuredThermalMesh
    series::ThermalSeries
    params::Params
    left_index::Int
    left::Array{Float64,3}
    right::Array{Float64,3}
    buffer::Array{Float64,3}
end

function ThermalFieldCache(mesh,series,p)
    length(series.times_s) >= 2 || throw(ArgumentError("Replay needs at least two thermal snapshots"))
    left = map_temperatures(mesh,series,p,first(series.times_s)*1000)
    right = map_temperatures(mesh,series,p,series.times_s[2]*1000)
    ThermalFieldCache(mesh,series,p,1,left,right,similar(left))
end

function cached_temperatures!(cache::ThermalFieldCache,time_ms)
    times = cache.series.times_s
    t = time_ms/1000
    first(times)-32eps(1.0) <= t <= last(times)+32eps(1.0) || throw(ArgumentError("Requested time lies outside available snapshots"))
    t = clamp(t,first(times),last(times))
    index = min(searchsortedlast(times,t),length(times)-1)
    if index != cache.left_index
        cache.left = index == cache.left_index+1 ? cache.right : map_temperatures(cache.mesh,cache.series,cache.params,times[index]*1000)
        cache.right = map_temperatures(cache.mesh,cache.series,cache.params,times[index+1]*1000)
        cache.left_index = index
    end
    fraction = (t-times[index])/(times[index+1]-times[index])
    Threads.@threads :static for index in eachindex(cache.buffer)
        cache.buffer[index] = (1-fraction)*cache.left[index] + fraction*cache.right[index]
    end
    cache.buffer
end

"""
    write_calculix_case(path, params; duration_ms=2.0, increment_ms=0.02)

Write a 2×2×4-element DC3D8 transient conduction specimen in SI units. Bottom
is held at 1500 K; top ramps from 1649 K to 1500 K; sides are insulated. Interior
nodal temperatures are solved. Density, heat capacity and conductivity are
illustrative constants, not a calibrated alloy. Returns mesh and initial nodal
snapshot because CalculiX need not print time zero.
"""
function write_calculix_case(path,p::Params; duration_ms=2.0,increment_ms=0.02)
    duration_ms > 0 && 0 < increment_ms <= duration_ms || throw(ArgumentError("Invalid solver times"))
    axes = [collect(range(0,L/1000;length=n)) for (L,n) in zip((p.lengthx,p.lengthy,p.lengthz),(3,3,5))]
    ids = reshape(collect(1:45),3,3,5)
    mesh = StructuredThermalMesh(axes...,ids)
    initial = [k==1 ? 1500.0 : 1649.0 for i in 1:3,j in 1:3,k in 1:5]
    duration,increment = duration_ms/1000,increment_ms/1000
    open(path,"w") do io
        println(io,"*HEADING\nAMMicrostructure one-way thermal integration smoke case (m, s, kg, K)\n*NODE,NSET=NALL")
        for k in 1:5,j in 1:3,i in 1:3
            println(io,join((ids[i,j,k],axes[1][i],axes[2][j],axes[3][k]),","))
        end
        println(io,"*ELEMENT,TYPE=DC3D8,ELSET=EALL")
        e=0
        for k in 1:4,j in 1:2,i in 1:2
            e+=1
            println(io,join((e,ids[i,j,k],ids[i+1,j,k],ids[i+1,j+1,k],ids[i,j+1,k],
                ids[i,j,k+1],ids[i+1,j,k+1],ids[i+1,j+1,k+1],ids[i,j+1,k+1]),","))
        end
        println(io,"*NSET,NSET=BOTTOM\n",join(vec(ids[:,:,1]),","))
        println(io,"*NSET,NSET=TOP\n",join(vec(ids[:,:,end]),","))
        println(io,"*MATERIAL,NAME=ILLUSTRATIVE\n*CONDUCTIVITY\n25.0\n*DENSITY\n7800.0\n*SPECIFIC HEAT\n600.0")
        println(io,"*SOLID SECTION,ELSET=EALL,MATERIAL=ILLUSTRATIVE\n*INITIAL CONDITIONS,TYPE=TEMPERATURE")
        for id in 1:45
            println(io,"$id,$(initial[id])")
        end
        println(io,"*AMPLITUDE,NAME=TOPCOOL,TIME=TOTAL TIME\n0.0,1649.0,$duration,1500.0")
        println(io,"*STEP,INC=10000\n*HEAT TRANSFER,DIRECT\n$increment,$duration")
        println(io,"*BOUNDARY\nBOTTOM,11,11,1500.0\n*BOUNDARY,AMPLITUDE=TOPCOOL\nTOP,11,11,1.0")
        println(io,"*NODE PRINT,NSET=NALL,FREQUENCY=1\nNT\n*NODE FILE,FREQUENCY=1\nNT\n*END STEP")
    end
    mesh,ThermalSeries([0.0],vec(ids),reshape(vec(initial),:,1))
end

"""Prescribed layer-top heating, followed by release (`nothing`), in ms/K."""
function layer_top_temperature(time_ms; pulse_ms=0.4,ramp_ms=0.2,base_K=1500.0,peak_K=1800.0)
    all(isfinite,(time_ms,pulse_ms,ramp_ms,base_K,peak_K)) && time_ms >= 0 &&
        0 < ramp_ms <= pulse_ms && 0 <= base_K < peak_K || throw(ArgumentError("Invalid layer heating schedule"))
    time_ms > pulse_ms && return nothing
    base_K + (peak_K-base_K)*min(time_ms/ramp_ms,1.0)
end

"""
Write a stock-ccx approximation of the supplied AscentAM test-deck pattern:
prescribe a hot layer top, release it, then cool by conduction and a top film.
This is a single existing solid layer, without activation or a moving laser.
"""
function write_calculix_layer_case(path,p::Params;duration_ms=2.0,increment_ms=0.02,cycles=1,localized=false)
    cycles in (1,2) || throw(ArgumentError("One or two heating cycles are supported"))
    period_ms = duration_ms/cycles
    period_ms > 0.4 && 0 < increment_ms <= 0.2 || throw(ArgumentError("Invalid solver times"))
    pulse_ms,ramp_ms = 0.4,0.2
    for interval in (pulse_ms,period_ms-pulse_ms)
        isapprox(interval/increment_ms,round(interval/increment_ms);atol=1e-8) ||
            throw(ArgumentError("Heating and cooling periods must contain whole increments"))
    end
    dims = localized ? (17,3,9) : (5,3,9)
    axes = [collect(range(0,L/1000;length=n)) for (L,n) in zip((p.lengthx,p.lengthy,p.lengthz),dims)]
    ids = reshape(collect(1:prod(dims)),dims...)
    mesh = StructuredThermalMesh(axes...,ids)
    initial = ThermalSeries([0.0],vec(ids),fill(1500.0,length(ids),1))
    open(path,"w") do io
        println(io,"*HEADING\nAscentAM-inspired layer-top heat/release/cool (m, s, kg, K); illustrative\n*NODE,NSET=NALL")
        for k in 1:dims[3],j in 1:dims[2],i in 1:dims[1]
            println(io,join((ids[i,j,k],axes[1][i],axes[2][j],axes[3][k]),","))
        end
        println(io,"*ELEMENT,TYPE=DC3D8,ELSET=EALL")
        e = 0
        top_elements = Int[]
        for k in 1:dims[3]-1,j in 1:dims[2]-1,i in 1:dims[1]-1
            e += 1
            println(io,join((e,ids[i,j,k],ids[i+1,j,k],ids[i+1,j+1,k],ids[i,j+1,k],
                ids[i,j,k+1],ids[i+1,j,k+1],ids[i+1,j+1,k+1],ids[i,j+1,k+1]),","))
            k == dims[3]-1 && push!(top_elements,e)
        end
        for (header,entries) in (("*ELSET,ELSET=ETOP",top_elements),
                                  ("*NSET,NSET=BOTTOM",vec(ids[:,:,1])),
                                  ("*NSET,NSET=TOP",vec(ids[:,:,end])))
            println(io,header)
            for block in Iterators.partition(entries,16)
                println(io,join(block,","))
            end
        end
        println(io,"*MATERIAL,NAME=ILLUSTRATIVE\n*CONDUCTIVITY\n25.0\n*DENSITY\n7800.0\n*SPECIFIC HEAT\n600.0")
        println(io,"*SOLID SECTION,ELSET=EALL,MATERIAL=ILLUSTRATIVE\n*INITIAL CONDITIONS,TYPE=TEMPERATURE\nNALL,1500.0")
        for cycle in 1:cycles
            start_ms = (cycle-1)*period_ms
            center_m = p.lengthx/1000*(cycle == 1 ? 0.375 : 0.625)
            for i in 1:dims[1]
                peak = localized ? 1500+300exp(-0.5*((axes[1][i]-center_m)/1e-5)^2) : 1800.0
                println(io,"*AMPLITUDE,NAME=HEAT$(cycle)X$(i),TIME=TOTAL TIME\n$(start_ms/1000),1500.0,$((start_ms+ramp_ms)/1000),$peak,$((start_ms+pulse_ms)/1000),$peak")
            end
        end
        for cycle in 1:cycles
            println(io,"*STEP,INC=10000\n*HEAT TRANSFER,DIRECT\n$(increment_ms/1000),$(pulse_ms/1000)")
            println(io,"*BOUNDARY,OP=NEW\nBOTTOM,11,11,1500.0")
            for i in 1:dims[1]
                println(io,"*BOUNDARY,AMPLITUDE=HEAT$(cycle)X$(i)")
                for j in 1:dims[2]
                    println(io,"$(ids[i,j,end]),11,11,1.0")
                end
            end
            # Remove the previous cooling film while the top is prescribed.
            cycle > 1 && println(io,"*FILM,OP=NEW\nETOP,F2,1500.0,0.0")
            println(io,"*NODE PRINT,NSET=NALL,FREQUENCY=1\nNT\n*NODE FILE,FREQUENCY=1\nNT\n*END STEP")
            println(io,"*STEP,INC=10000\n*HEAT TRANSFER,DIRECT\n$(increment_ms/1000),$((period_ms-pulse_ms)/1000)")
            println(io,"*BOUNDARY,OP=NEW\nBOTTOM,11,11,1500.0\n*FILM,OP=NEW\nETOP,F2,1500.0,100000.0")
            println(io,"*NODE PRINT,NSET=NALL,FREQUENCY=1\nNT\n*NODE FILE,FREQUENCY=1\nNT\n*END STEP")
        end
    end
    mesh,initial
end

function write_slice_frame(io, state, temp, y_index)
    p = state.params
    for k in 1:p.nz, i in 1:p.nx
        @printf(io,"%d,%.10g,%d,%d,%d,%d,%d,%.10g\n",state.step,state.time_ms,
            i,y_index,k,state.ids[i,y_index,k],state.states[i,y_index,k],temp[i,y_index,k])
    end
end

function evolve_thermal_history(mesh,series,c,output_dir;solver_executed=false,source_kind="unknown")
    p = parameters(c)
    output = abspath(output_dir)
    mkpath(output)
    open(io->TOML.print(io,c;sorted=true),joinpath(output,"config.toml"),"w")
    initialization = get(c,"ca_initialization","bottom-layer")
    thermal_cache = ThermalFieldCache(mesh,series,p)
    state = SimulationState(p,cached_temperatures!(thermal_cache,0.0);seed=c["seed"],initialization,
        max_cells=get(c,"max_cells",250_000))
    rows = [(step=0,time_ms=0.0,assigned_fraction=count(>(0),state.ids)/length(state.ids),grains=state.num_grains,captures=0)]
    melted_cells = falses(size(state.ids))
    max_unassigned = count(iszero,state.ids)
    write_vtk(joinpath(output,"initial.vtk"),p,state.ids,state.states,state.previous_temperature)
    snapshots = [(time=0.0,fraction=first(rows).assigned_fraction,ids=copy(state.ids))]
    # Stream physical x-z slices at the configured cadence; retain no full 3D history.
    slice_y = cld(p.ny,2)
    slice_stride = get(c,"slice_stride",1)
    slice_stride isa Integer && slice_stride > 0 || throw(ArgumentError("slice_stride must be a positive integer"))
    slice_frames = 1
    cold_skipped_steps = 0
    open(joinpath(output,"slice_history.csv"),"w") do slice_io
      println(slice_io,"step,time_ms,x_index,y_index,z_index,grain_id,cell_state,temperature_K")
      write_slice_frame(slice_io,state,state.previous_temperature,slice_y)
      open(joinpath(output,"kernel.log"),"w") do io
        redirect_stdout(io) do
            for n in 1:c["steps"]
                temp = cached_temperatures!(thermal_cache,n*p.dt)
                melted_cells .|= temp .>= p.Tliquidus
                logvars = step!(state,temp;skip_cold=get(c,"skip_cold",true))
                cold_skipped_steps += logvars["cold_shortcut"]
                if n % slice_stride == 0 || n == c["steps"]
                    write_slice_frame(slice_io,state,temp,slice_y)
                    slice_frames += 1
                end
                if get(c,"thermal_case","") == "two-pass-heating"
                    n == round(Int,2.0/p.dt) && write_vtk(joinpath(output,"first_cycle.vtk"),p,state.ids,state.states,temp)
                    n == round(Int,2.4/p.dt) && write_vtk(joinpath(output,"second_peak.vtk"),p,state.ids,state.states,temp)
                end
                max_unassigned = max(max_unassigned,count(iszero,state.ids))
                push!(rows,(step=n,time_ms=state.time_ms,assigned_fraction=count(>(0),state.ids)/length(state.ids),grains=state.num_grains,captures=Int(logvars["numOctahedralGrowth"])))
                if n == get(c,"snapshot_step",cld(c["steps"],2)) || n == c["steps"]
                    push!(snapshots,(time=state.time_ms,fraction=last(rows).assigned_fraction,ids=copy(state.ids)))
                    n < c["steps"] && write_vtk(joinpath(output,"peak.vtk"),p,state.ids,state.states,temp)
                end
            end
        end
      end
    end
    write_vtk(joinpath(output,"final.vtk"),p,state.ids,state.states,state.previous_temperature)
    open(joinpath(output,"history.csv"),"w") do io
        println(io,"step,time_ms,assigned_fraction,grains,octahedral_captures")
        for row in rows
            println(io,join(values(row),","))
        end
    end
    preview_rows = [(time=r.time_ms,fraction=r.assigned_fraction) for r in rows]
    write_preview(joinpath(output,"preview.svg"),p,snapshots,preview_rows,c["seed"];thermal_label="Saved FEM temperatures")
    summary = Dict("solver_executed"=>solver_executed,"thermal_source"=>source_kind,
        "solver_snapshots"=>length(series.times_s)-1,"cells"=>length(state.ids),"steps"=>state.step,
        "time_ms"=>state.time_ms,"octahedral_captures"=>sum(r.captures for r in rows),
        "initial_grains"=>p.numGrains,"final_grains"=>state.num_grains,
        "nx"=>p.nx,"ny"=>p.ny,"nz"=>p.nz,"dx_um"=>p.dx*1000,"dt_ms"=>p.dt,
        "resolution"=>get(c,"resolution","quick"),"slice_frames"=>slice_frames,"slice_stride"=>slice_stride,
        "worker_threads"=>Threads.nthreads(),"cold_skipped_steps"=>cold_skipped_steps,
        "active_steps"=>state.step-cold_skipped_steps,"rng_contract"=>"Xoshiro seed and cell-index capture streams; worker-count independent",
        "growth_rule"=>"active solid with a mushy cell in its 26 first-Moore neighbours",
        "cold_skip_policy"=>"both endpoint fields below solidus and all cells assigned active solid: freeze octahedra, skip microsteps, advance temperature/time without RNG draws",
        "final_slice_grains"=>length(filter(>(0),unique(state.ids[:,slice_y,:]))),
        "assigned_fraction"=>last(rows).assigned_fraction,"finite_octahedra"=>all(isfinite,state.octa),
        "temperature_min_K"=>minimum(series.temperature_K),"temperature_max_K"=>maximum(series.temperature_K),
        "ca_initialization"=>initialization,"melted_cells"=>count(melted_cells),"max_unassigned_cells"=>max_unassigned,
        "snapshot_time_ms"=>get(c,"snapshot_step",cld(c["steps"],2))*p.dt,
        "slice_y_index"=>slice_y,"slice_y_um"=>(slice_y-0.5)*p.dx*1000)
    for key in ("heating_intervals_ms","thermal_case","thermal_dataset")
        haskey(c,key) && (summary[key] = c[key])
    end
    open(io->TOML.print(io,summary;sorted=true),joinpath(output,"summary.toml"),"w")
    println("Nodal thermal history → trilinear cell temperatures → original CA complete: ",output)
    (;state,series,mesh,summary,rows)
end

file_sha256(path) = bytes2hex(open(sha256,path))
const THERMAL_DATA_FILES = ("thermal.dat","thermal.inp","mesh.toml","initial.toml","config.toml")

"""Write portable mesh/configuration sidecars and hashes beside an existing raw solver .dat/.inp."""
function write_thermal_dataset(directory,mesh,initial,c;provenance)
    mesh_data = Dict("coordinate_unit"=>"m","x_m"=>mesh.x_m,"y_m"=>mesh.y_m,"z_m"=>mesh.z_m,
        "node_ids"=>vec(mesh.node_ids),"node_order"=>"x-fastest")
    initial_data = Dict("temperature_unit"=>"K","time_s"=>0.0,"node_ids"=>initial.node_ids,
        "temperature_K"=>vec(initial.temperature_K),"ca_initialization"=>get(c,"ca_initialization","bottom-layer"),
        "ca_seed"=>c["seed"])
    for (name,data) in (("mesh.toml",mesh_data),("initial.toml",initial_data),("config.toml",c))
        open(io->TOML.print(io,data;sorted=true),joinpath(directory,name),"w")
    end
    metadata = merge(Dict("schema_version"=>1,"julia_version"=>string(VERSION)),provenance)
    metadata["sha256"] = Dict(name=>file_sha256(joinpath(directory,name)) for name in THERMAL_DATA_FILES)
    open(io->TOML.print(io,metadata;sorted=true),joinpath(directory,"provenance.toml"),"w")
    metadata
end

function load_thermal_dataset(directory)
    metadata = TOML.parsefile(joinpath(directory,"provenance.toml"))
    metadata["schema_version"] == 1 || error("Unsupported thermal dataset schema")
    for name in THERMAL_DATA_FILES
        file_sha256(joinpath(directory,name)) == metadata["sha256"][name] || error("Thermal dataset checksum mismatch: $name")
    end
    c = load_demo_config(joinpath(directory,"config.toml"))
    m = TOML.parsefile(joinpath(directory,"mesh.toml"))
    m["coordinate_unit"] == "m" && m["node_order"] == "x-fastest" || error("Unsupported mesh units/layout")
    axes = (m["x_m"],m["y_m"],m["z_m"])
    mesh = StructuredThermalMesh(axes...,reshape(m["node_ids"],length.(axes)...))
    initial = TOML.parsefile(joinpath(directory,"initial.toml"))
    initial["time_s"] == 0 && initial["temperature_unit"] == "K" || error("Expected initial snapshot at 0 s in K")
    initial["ca_seed"] == c["seed"] || error("Initial substrate seed differs from configuration")
    solved = read_calculix_dat(joinpath(directory,"thermal.dat"))
    solved.node_ids == initial["node_ids"] || error("Solver node IDs differ from initial snapshot")
    first(solved.times_s) > 0 || error("Expected strictly positive solver output times")
    series = ThermalSeries(vcat(0.0,solved.times_s),solved.node_ids,hcat(initial["temperature_K"],solved.temperature_K))
    (;mesh,series,config=c,provenance=metadata)
end

"""Apply a named CA sampling preset without modifying the archived solver dataset."""
function thermal_replay_config(config,resolution)
    profiles = TOML.parsefile(joinpath(@__DIR__,"..","examples","fem_replay","resolutions.toml"))
    haskey(profiles,resolution) || throw(ArgumentError("Resolution must be showcase, fine or quick"))
    profile = profiles[resolution]
    c = deepcopy(config)
    spatial,temporal = profile["spatial_refinement"],profile["temporal_refinement"]
    for axis in ("nx","ny","nz")
        c["grid"][axis] *= spatial
    end
    c["grid"]["dx_mm"] /= spatial
    c["grid"]["substrate_grains"] = get(profile,"substrate_grains",c["grid"]["substrate_grains"])
    c["dt_ms"] /= temporal
    c["steps"] *= temporal
    haskey(c,"snapshot_step") && (c["snapshot_step"] *= temporal)
    c["resolution"] = resolution
    c["slice_stride"] = get(profile,"slice_stride",1)
    haskey(profile,"max_cells") && (c["max_cells"] = profile["max_cells"])
    parameters(c) # Check the allocation guard before any large arrays are built.
    c
end

"""
    run_thermal_replay(; dataset_dir="data/two_pass_heating_demo", output_dir="output/fem_replay", resolution="showcase")

Run the original CA from a portable nodal thermal dataset without launching or
requiring CalculiX. Validates file SHA256 checksums, units and mesh/time coverage.
The saved seed/configuration reconstruct the substrate initialization; exact
random-stream reproducibility is supported on the same Julia version.
The default showcase preset samples the bundled domain with 6,635,520 cells and
64 substrate grains. The fine preset uses 245,760 cells and 32 grains; quick
retains the small archived CA settings. Showcase and fine retain every fourth
CA step as a slice and every step in the scalar history.
"""
function run_thermal_replay(;dataset_dir=joinpath(@__DIR__,"..","data","two_pass_heating_demo"),
                            output_dir=joinpath(@__DIR__,"..","output","fem_replay"),resolution="showcase")
    dataset = load_thermal_dataset(dataset_dir)
    config = thermal_replay_config(dataset.config,resolution)
    config["thermal_dataset"] = relpath(abspath(dataset_dir),joinpath(@__DIR__,".."))
    result = evolve_thermal_history(dataset.mesh,dataset.series,config,output_dir;
        source_kind=dataset.provenance["source_kind"])
    println("Replayed saved temperatures without a solver: ",abspath(dataset_dir))
    result
end

"""Run stock ccx, save a reusable nodal thermal dataset, then advance the original CA."""
function run_calculix_demo(;ccx=get(ENV,"CCX","ccx"),output_dir=joinpath(@__DIR__,"..","output","calculix"),case="two-pass-heating",resolution="showcase")
    resolution in ("showcase","fine","quick") || throw(ArgumentError("Resolution must be showcase, fine or quick"))
    executable = Sys.which(ccx)
    executable === nothing && error("CalculiX executable not found: $ccx. Install ccx or set CCX=/absolute/path/to/ccx; see docs/fem-integration.md")
    c = load_demo_config(joinpath(@__DIR__,"..","examples","directional_solidification","config.toml"))
    case in ("conduction","layer-heating","two-pass-heating") || throw(ArgumentError("Unknown thermal case: $case"))
    c["grid"]["nx"],c["grid"]["ny"],c["grid"]["nz"] = case == "conduction" ? (12,4,12) : (32,6,20)
    case != "conduction" && (c["ca_initialization"] = "solid-substrate")
    case != "conduction" && (c["snapshot_step"] = 20)
    c["thermal_case"] = case
    c["heating_intervals_ms"] = case == "two-pass-heating" ? [[0.0,0.4],[2.0,2.4]] : case == "layer-heating" ? [[0.0,0.4]] : Vector{Float64}[]
    case == "two-pass-heating" && (c["steps"] = 200)
    p = parameters(c)
    output = abspath(output_dir)
    mkpath(output)
    work = mktempdir(output;prefix="solver-",cleanup=false)
    writer = case == "conduction" ? write_calculix_case : write_calculix_layer_case
    kwargs = case == "two-pass-heating" ? (;cycles=2,localized=true) : (;)
    mesh,initial = writer(joinpath(work,"thermal.inp"),p;duration_ms=c["steps"]*p.dt,increment_ms=p.dt,kwargs...)
    open(joinpath(work,"solver.log"),"w") do io
        command = Cmd(Cmd([executable,"-i","thermal"]);dir=work)
        try
            run(pipeline(command,stdout=io,stderr=io))
        catch exception
            exception isa ProcessFailedException || rethrow()
            error("CalculiX failed before coupling. Inspect $(joinpath(work, "solver.log")) for solver or missing shared-library errors; no new coupled summary was written.")
        end
    end
    solver_log = read(joinpath(work,"solver.log"),String)
    version = match(r"Version\s+([0-9]+(?:\.[0-9]+)+)",solver_log)
    provenance = Dict("source_kind"=>"calculix","solver"=>"CalculiX CrunchiX",
        "solver_version"=>version === nothing ? "unknown (see solver.log)" : version[1],
        "solver_binary_sha256"=>file_sha256(executable),"generated_unix_s"=>time(),"solver_command"=>"ccx -i thermal",
        "thermal_case"=>case,
        "description"=>case == "conduction" ? "Actual transient conduction output; 16 DC3D8 elements, 45 nodes; illustrative material" :
            case == "two-pass-heating" ? "Two actual localized top-heating/release/cooling cycles; 256 DC3D8 elements, 459 nodes; Gaussian top temperature centered at 24 and 40 micrometres, sigma 10 micrometres; prescribed stationary pulses, no moving laser" :
            "Actual layer-top heating/release/cooling; 64 DC3D8 elements, 135 nodes; AscentAM test-deck pattern with illustrative values")
    write_thermal_dataset(work,mesh,initial,c;provenance)
    dataset = load_thermal_dataset(work)
    config = thermal_replay_config(dataset.config,resolution)
    result = evolve_thermal_history(dataset.mesh,dataset.series,config,output;
        solver_executed=true,source_kind="calculix")
    result.summary["solver_executable"] = executable
    result.summary["solver_directory"] = work
    result.summary["solver_version"] = provenance["solver_version"]
    open(io->TOML.print(io,result.summary;sorted=true),joinpath(output,"summary.toml"),"w")
    println("Saved reusable solver dataset: ",work)
    merge(result,(dataset_dir=work,))
end
