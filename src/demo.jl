
using Random, TOML, Printf

function load_demo_config(path)
    c = TOML.parsefile(path)
    g, t, m = c["grid"], c["thermal"], c["material"]
    all(g[k] isa Integer && g[k] > 0 for k in ("nx", "ny", "nz", "substrate_grains")) ||
        error("Grid dimensions and substrate_grains must be positive integers")
    check_grid(g["nx"],g["ny"],g["nz"],g["substrate_grains"];max_cells=get(c,"max_cells",250_000))
    c["steps"] > 0 && c["dt_ms"] > 0 && g["dx_mm"] > 0 || error("steps, dt_ms and dx_mm must be positive")
    m["liquidus_K"] > m["solidus_K"] + 1 || error("Liquidus must exceed solidus by more than 1 K")
    m["nucleation_sigma_K"] > 0 && m["nucleation_N0"] >= 0 || error("Invalid nucleation parameters")
    t["cooling_rate_K_per_ms"] > 0 || error("Cooling rate must be positive")
    all(isfinite, [v for d in (g,t,m) for v in values(d)]) || error("Configuration values must be finite")
    return c
end

function parameters(c)
    g, m = c["grid"], c["material"]
    nx, ny, nz, dx = g["nx"], g["ny"], g["nz"], g["dx_mm"]
    check_grid(nx,ny,nz,g["substrate_grains"];max_cells=get(c,"max_cells",250_000))
    return Params(Tliquidus=m["liquidus_K"], Tsolidus=m["solidus_K"], tau=1.0,
        lengthx=nx*dx, lengthy=ny*dx, lengthz=nz*dx, nx=nx, ny=ny, nz=nz,
        resolution=max(nx,ny,nz), dx=dx, dt=c["dt_ms"], numGrains=g["substrate_grains"],
        minGrainSize=dx, D_L=m["D_L"], Gamma=m["Gamma"], m_L=m["m_L"],
        k_0=m["k_0"], C_0=m["C_0"], dTsigma=m["nucleation_sigma_K"],
        dTc=m["nucleation_mean_K"], N0=m["nucleation_N0"],
        configName="examples/directional_solidification/config.toml", filename="standalone", multilayer=false)
end

function temperatures(p, c, time_ms)
    t = c["thermal"]
    values = [p.Tliquidus + t["superheat_K"] + t["gradient_K_per_mm"]*(k-0.5)*p.dx -
              t["cooling_rate_K_per_ms"]*time_ms for i in 1:p.nx, j in 1:p.ny, k in 1:p.nz]
    # A prescribed cold bottom layer supplies epitaxial substrate grains.
    values[:,:,1] .= p.Tsolidus - 5
    return values
end

function initialize(p)
    ids = zeros(Int, p.nx, p.ny, p.nz)
    states = fill(3, size(ids))
    octa = zeros(4, size(ids)...)
    grains = zeros(3, length(ids) + p.numGrains)
    grains[:,1:p.numGrains] .= 90 .* rand(3,p.numGrains) .- 45
    sites = shuffle(collect(CartesianIndices((p.nx,p.ny))))[1:p.numGrains]
    for j in 1:p.ny, i in 1:p.nx
        ids[i,j,1] = argmin([sum(abs2, Tuple(CartesianIndex(i,j)-s)) for s in sites])
        states[i,j,1] = 1
    end
    return ids, states, octa, grains
end

function write_vtk(path, p, ids, states, temp)
    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0\nStandalone thesis CA demonstration\nASCII\nDATASET STRUCTURED_POINTS")
        println(io, "DIMENSIONS $(p.nx+1) $(p.ny+1) $(p.nz+1)\nORIGIN 0 0 0\nSPACING $(p.dx) $(p.dx) $(p.dx)")
        println(io, "CELL_DATA $(length(ids))")
        for (name, kind, data) in (("grain_id", "int", ids), ("cell_state", "int", states), ("temperature_K", "double", temp))
            println(io, "SCALARS $name $kind 1\nLOOKUP_TABLE default")
            foreach(x -> println(io,x), data)
        end
    end
end

grain_color(id) = id == 0 ? "#101c30" : "hsl($(mod(id*137.508,360)),64%,61%)"

function write_preview(path, p, snapshots, rows, seed;thermal_label="Prescribed cooling")
    open(path,"w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="610" viewBox="0 0 1100 610">
        <rect width="1100" height="610" fill="#08121f"/>
        <g font-family="Arial,sans-serif" fill="#edf4fc">
        <text x="34" y="42" font-size="25" font-weight="bold">Directional solidification · original Julia cellular automaton</text>
        <text x="34" y="70" font-size="15" fill="#afc4d9">$thermal_label · decentered octahedral capture · seed $seed · illustrative, uncalibrated</text>""")
        scale = min(310/p.nx, 270/p.nz)
        for (panel,snap) in enumerate(snapshots)
            x0 = 34 + (panel-1)*358
            y0 = 118
            @printf(io,"<text x=\"%d\" y=\"103\" font-size=\"16\">t = %.2f ms · %.1f%% assigned</text>\n",x0,snap.time,100snap.fraction)
            for k in 1:p.nz, i in 1:p.nx
                id = snap.ids[i,cld(p.ny,2),k]
                @printf(io,"<rect x=\"%.3f\" y=\"%.3f\" width=\"%.3f\" height=\"%.3f\" fill=\"%s\"/>\n",x0+(i-1)*scale,y0+(p.nz-k)*scale,scale+0.02,scale+0.02,grain_color(id))
            end
            @printf(io,"<text x=\"%d\" y=\"390\" font-size=\"13\" fill=\"#afc4d9\">x–z slice at y = %.3f mm · width %.3f mm</text>\n", x0,(cld(p.ny,2)-0.5)*p.dx,p.lengthx)
        end
        println(io,"<text x=\"34\" y=\"427\" font-size=\"16\">Assigned grain fraction (whole 3D domain)</text>")
        println(io,"<path d=\"M 65 451 L 65 550 L 1045 550\" fill=\"none\" stroke=\"#67819e\"/>")
        points = join([@sprintf("%.2f,%.2f",65+980*r.time/last(rows).time,550-99*r.fraction) for r in rows]," ")
        println(io,"<polyline points=\"$points\" fill=\"none\" stroke=\"#69d4d2\" stroke-width=\"3\"/>")
        println(io,"<text x=\"26\" y=\"457\" font-size=\"12\">100%</text><text x=\"44\" y=\"554\" font-size=\"12\">0</text>")
        @printf(io,"<text x=\"65\" y=\"573\" font-size=\"12\">0 ms</text><text x=\"990\" y=\"573\" font-size=\"12\">%.2f ms</text>\n",last(rows).time)
        println(io,"<text x=\"34\" y=\"600\" font-size=\"12\" fill=\"#afc4d9\">Colors identify grains; dark cells have no grain yet. Bottom layer is prescribed substrate. Assignment is a CA state, not thermodynamic solid fraction.</text></g></svg>")
    end
end

function run_demo(; config_path=joinpath(@__DIR__,"..","examples","directional_solidification","config.toml"), output_dir=joinpath(@__DIR__,"..","output","demo"))
    c = load_demo_config(config_path)
    p = parameters(c)
    simulation = SimulationState(p,temperatures(p,c,0.0);seed=c["seed"],max_cells=get(c,"max_cells",250_000))
    ids, states, octa, grains = simulation.ids, simulation.states, simulation.octa, simulation.grains
    mkpath(output_dir)
    open(joinpath(output_dir,"config.toml"),"w") do io
        TOML.print(io,c; sorted=true)
    end
    previous = temperatures(p,c,0.0)
    ngrains = p.numGrains
    rows = [(step=0,time=0.0,fraction=count(>(0),ids)/length(ids),grains=ngrains,nucleated=0,captured=0)]
    snapshots = [(time=0.0,fraction=rows[1].fraction,ids=copy(ids))]
    write_vtk(joinpath(output_dir,"initial.vtk"),p,ids,states,previous)
    println("Running $(p.nx) × $(p.ny) × $(p.nz) cells, $(c["steps"]) steps, seed $(c["seed"])…")
    open(joinpath(output_dir,"kernel.log"),"w") do logio
        redirect_stdout(logio) do
            for step in 1:c["steps"]
                temp = temperatures(p,c,step*p.dt)
                logvars = step!(simulation,temp)
                ngrains = simulation.num_grains
                push!(rows,(step=step,time=step*p.dt,fraction=count(>(0),ids)/length(ids),grains=ngrains,
                    nucleated=Int(logvars["numCellsNucleated"]),captured=Int(logvars["numOctahedralGrowth"])))
                if step == cld(c["steps"],2) || step == c["steps"]
                    push!(snapshots,(time=step*p.dt,fraction=last(rows).fraction,ids=copy(ids)))
                    write_vtk(joinpath(output_dir,step==c["steps"] ? "final.vtk" : "midpoint.vtk"),p,ids,states,temp)
                end
                previous = temp
            end
        end
    end
    open(joinpath(output_dir,"history.csv"),"w") do io
        println(io,"step,time_ms,assigned_fraction,grains,nucleated,octahedral_captures")
        for r in rows
            println(io,join(values(r),","))
        end
    end
    open(joinpath(output_dir,"grains.csv"),"w") do io
        println(io,"grain_id,phi1_deg,phi2_deg,phi3_deg,cell_count")
        for id in 1:ngrains
            println(io,join((id,grains[:,id]...,count(==(id),ids)),","))
        end
    end
    write_preview(joinpath(output_dir,"preview.svg"),p,snapshots,rows,c["seed"])
    summary = Dict("seed"=>c["seed"],"cells"=>length(ids),"steps"=>c["steps"],"julia_version"=>string(VERSION),
        "time_ms"=>last(rows).time,"initial_grains"=>p.numGrains,"final_grains"=>ngrains,
        "assigned_fraction"=>last(rows).fraction,"nucleated_grains"=>sum(r.nucleated for r in rows),
        "octahedral_captures"=>sum(r.captured for r in rows),"threads"=>Threads.nthreads(),
        "finite_octahedra"=>all(isfinite,octa),"valid_grain_ids"=>all(id->0<=id<=ngrains,ids))
    open(joinpath(output_dir,"summary.toml"),"w") do io
        TOML.print(io,summary; sorted=true)
    end
    @printf("Done: %.1f%% assigned, %d grains, %d octahedral captures.\n",100last(rows).fraction,ngrains,summary["octahedral_captures"])
    println("Outputs: ",abspath(output_dir))
    return (;summary,ids,states,octa,grains,rows)
end


