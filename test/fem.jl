@testset "CalculiX text output contract" begin
    series = read_calculix_dat(joinpath(@__DIR__,"fixtures","temperatures.dat"))
    @test series.times_s == [0.001,0.002]
    @test series.node_ids == [10,30]
    @test series.temperature_K == [1600.0 1500.0;1620.0 1520.0]
    @test CA.nodal_temperatures(series,0.0015) ≈ [1550,1570]
    @test CA.nodal_temperatures(series,0.002) == [1500,1520]
    @test_throws ArgumentError CA.nodal_temperatures(series,0.0)
    @test_throws ArgumentError ThermalSeries([0.0,0.0],[1],ones(1,2))
    @test_throws ArgumentError ThermalSeries([0.0],[1],fill(NaN,1,1))
    @test_throws ArgumentError ThermalSeries([0.0],[1,1],ones(2,1))
    mktempdir() do dir
        path = joinpath(dir,"bad.dat")
        write(path,"temperatures for set NALL and time 1.0\n 1 1600\n 1 1500\n")
        @test_throws ErrorException read_calculix_dat(path)
        write(path,"temperatures for set NALL and time 1.0\n 1 1600\ntemperatures for set NALL and time 2.0\n 2 1500\n")
        @test_throws ErrorException read_calculix_dat(path)
        write(path,"temperatures for set NALL and time 1.0\n 1 1600\ntemperatures for set NALL and time 0.5\n 1 1500\n")
        @test_throws ArgumentError read_calculix_dat(path)
    end
end

@testset "Spatial and temporal mapping, units and bounds" begin
    c = CA.load_demo_config(CONFIG)
    c["grid"]["nx"],c["grid"]["ny"],c["grid"]["nz"] = 4,3,5
    p = CA.parameters(c)
    axes = [[0.0,L/2000,L/1000] for L in (p.lengthx,p.lengthy,p.lengthz)]
    ids = reshape(collect(101:127),3,3,3)
    mesh = StructuredThermalMesh(axes...,ids)
    # A multilinear field checks all eight shape functions, SI conversion,
    # interpolation in time, and lookup by arbitrary unsorted node IDs.
    field(x,y,z,t) = 1500 + 10x/p.lengthx + 20y/p.lengthy + 30z/p.lengthz +
        40x*y/(p.lengthx*p.lengthy) + 50x*y*z/(p.lengthx*p.lengthy*p.lengthz) - 10t
    nodal(t) = vec([field(1000x,1000y,1000z,t) for x in axes[1],y in axes[2],z in axes[3]])
    order = reverse(1:27)
    series = ThermalSeries([0,0.002],vec(ids)[order],hcat(nodal(0),nodal(2))[order,:])
    mapped = map_temperatures(mesh,series,p,0.7)
    expected = [field((i-0.5)*p.dx,(j-0.5)*p.dx,(k-0.5)*p.dx,0.7) for i in 1:p.nx,j in 1:p.ny,k in 1:p.nz]
    @test mapped ≈ expected atol=1e-10
    cache = CA.ThermalFieldCache(mesh,series,p)
    for time_ms in (0.0,0.7,2.0,0.3)
        @test CA.cached_temperatures!(cache,time_ms) ≈ map_temperatures(mesh,series,p,time_ms) atol=1e-10
    end
    @test_throws ArgumentError CA.cached_temperatures!(cache,2.1)
    three = ThermalSeries([0,0.001,0.002],series.node_ids,hcat(nodal(0),nodal(1),nodal(2))[order,:])
    cache3 = CA.ThermalFieldCache(mesh,three,p)
    for time_ms in (0.4,1.0,1.6,2.0,0.2)
        @test CA.cached_temperatures!(cache3,time_ms) ≈ map_temperatures(mesh,three,p,time_ms) atol=1e-10
    end
    @test_throws ArgumentError map_temperatures(mesh,series,p,2.1)
    missing = ThermalSeries([0,0.002],vec(ids)[1:26],hcat(nodal(0),nodal(2))[1:26,:])
    @test_throws ArgumentError map_temperatures(mesh,missing,p,1.0)
    @test_throws ArgumentError StructuredThermalMesh([0.0,0.0],axes[2],axes[3],ids)
    tiny = StructuredThermalMesh([0.0,1e-9],axes[2],axes[3],ids[1:2,:,:])
    @test_throws ArgumentError map_temperatures(tiny,series,p,1.0)
end

@testset "Generated deck and mapped thermal field drive original CA" begin
    mktempdir() do dir
        c = CA.load_demo_config(CONFIG)
        c["grid"]["nx"],c["grid"]["ny"],c["grid"]["nz"] = 8,3,8
        c["grid"]["substrate_grains"] = 4
        p = CA.parameters(c)
        deck = joinpath(dir,"thermal.inp")
        mesh,initial = write_calculix_case(deck,p)
        @test length(mesh.node_ids) == 45
        @test occursin("*ELEMENT,TYPE=DC3D8",read(deck,String))
        @test occursin("*HEAT TRANSFER,DIRECT",read(deck,String))
        @test initial.temperature_K[1] == 1500
        @test initial.temperature_K[end] == 1649
        # Explicitly synthetic history to test the coupling path without ccx.
        series = ThermalSeries([0.0,0.002],initial.node_ids,hcat(initial.temperature_K,fill(1500.0,45)))
        firsttemp = map_temperatures(mesh,series,p,0.0)
        state = SimulationState(p,firsttemp;seed=2026)
        firsttemp .= 999.0
        @test state.previous_temperature[1] != 999
        @test_throws DimensionMismatch step!(state,zeros(1,1,1))
        @test_throws ArgumentError step!(state,fill(NaN,size(state.ids)))
        @test_throws ArgumentError step!(state,state.previous_temperature;dt_ms=0)
        captures=0
        open(joinpath(dir,"kernel.log"),"w") do io
            redirect_stdout(io) do
                for n in 1:100
                    logvars = step!(state,map_temperatures(mesh,series,p,n*p.dt))
                    captures += Int(logvars["numOctahedralGrowth"])
                end
            end
        end
        @test state.step == 100
        @test state.time_ms ≈ 2.0
        @test captures > 0
        @test count(>(0),state.ids) > p.nx*p.ny
        @test all(isfinite,state.octa)
        @test all(id->0 <= id <= state.num_grains,state.ids)
        @test state.previous_temperature ≈ fill(1500,size(state.ids))
    end
end

@testset "Layer heating schedule and solid substrate initialization" begin
    @test CA.layer_top_temperature(0.0) == 1500.0
    @test CA.layer_top_temperature(0.1) == 1650.0
    @test CA.layer_top_temperature(0.2) == CA.layer_top_temperature(0.4) == 1800.0
    @test CA.layer_top_temperature(0.4001) === nothing # Top temperature constraint is released.
    @test_throws ArgumentError CA.layer_top_temperature(-0.1)
    @test_throws ArgumentError CA.layer_top_temperature(0.1;ramp_ms=0.5)
    c = CA.load_demo_config(CONFIG)
    c["grid"]["nx"],c["grid"]["ny"],c["grid"]["nz"] = 8,3,8
    p = CA.parameters(c)
    solid = SimulationState(p,fill(1500.0,p.nx,p.ny,p.nz);initialization="solid-substrate")
    @test all(>(0),solid.ids) && all(==(1),solid.states)
    @test solid.ids[:,:,end] == solid.ids[:,:,1]
    @test length(unique(solid.ids)) == p.numGrains
    @test_throws ArgumentError SimulationState(p,fill(Float64(p.Tsolidus),size(solid.ids));initialization="solid-substrate")
    @test_throws ArgumentError SimulationState(p,fill(1500.0,size(solid.ids));initialization="unknown")
    mktempdir() do dir
        mesh,initial = CA.write_calculix_layer_case(joinpath(dir,"thermal.inp"),p)
        @test length(mesh.node_ids) == 135
        @test all(==(1500.0),initial.temperature_K)
        @test_throws ArgumentError CA.write_calculix_layer_case(joinpath(dir,"bad.inp"),p;increment_ms=0.03)
        mesh2,initial2 = CA.write_calculix_layer_case(joinpath(dir,"two.inp"),p;duration_ms=4.0,cycles=2,localized=true)
        @test length(mesh2.node_ids) == 459
        @test all(==(1500.0),initial2.temperature_K)
        deck = read(joinpath(dir,"two.inp"),String)
        @test count(==("*STEP,INC=10000"),split(deck,'\n')) == 4
        @test occursin("*AMPLITUDE,NAME=HEAT2X1,TIME=TOTAL TIME",deck)
        @test all(length(split(line,',')) <= 16 for line in split(deck,'\n'))
        @test_throws ArgumentError CA.write_calculix_layer_case(joinpath(dir,"bad.inp"),p;cycles=0)
    end
end

@testset "Cold shortcut preserves geometry, RNG and later remelting" begin
    c = CA.load_demo_config(CONFIG)
    c["grid"]["nx"],c["grid"]["ny"],c["grid"]["nz"] = 8,3,8
    c["grid"]["substrate_grains"] = 4
    p = CA.parameters(c)
    initial = fill(1500.0,p.nx,p.ny,p.nz)
    fast = SimulationState(p,initial;initialization="solid-substrate",seed=2026)
    reference = deepcopy(fast)
    mktemp() do _,io
      redirect_stdout(io) do
        for offset in (0.0,10.0,20.0,5.0)
            temp = initial .+ offset
            log = step!(fast,temp)
            step!(reference,temp;skip_cold=false)
            @test log["cold_shortcut"]
            @test log["substeps"] == 0
            @test all(iszero,fast.octa)
            @test fast.ids == reference.ids
            @test fast.states == reference.states
            @test fast.octa == reference.octa
            @test fast.previous_temperature == reference.previous_temperature
            @test fast.time_ms == reference.time_ms
            @test rand(copy(fast.rng),UInt64) == rand(copy(reference.rng),UInt64)
        end
        # A future remelt must see identical frozen cold octahedron lengths.
        for top_temperature in (1660.0,1580.0,1540.0,1500.0,1500.0)
            temp = copy(initial)
            temp[:,:,5:end] .= top_temperature
            step!(fast,temp)
            step!(reference,temp;skip_cold=false)
            @test fast.ids == reference.ids
            @test fast.states == reference.states
            @test fast.octa == reference.octa
            @test fast.grains == reference.grains
            @test fast.num_grains == reference.num_grains
            @test rand(copy(fast.rng),UInt64) == rand(copy(reference.rng),UInt64)
        end
        incomplete = SimulationState(p,initial;initialization="solid-substrate")
        incomplete.ids[4,2,4],incomplete.states[4,2,4] = 0,2
        @test !step!(incomplete,initial)["cold_shortcut"]
        warm = SimulationState(p,initial;initialization="solid-substrate")
        warm.previous_temperature[4,2,4] = p.Tsolidus
        @test !step!(warm,initial)["cold_shortcut"]
      end
    end
end

@testset "Octahedra grow only beside first-Moore mushy neighbours" begin
    c = CA.load_demo_config(CONFIG)
    c["grid"]["nx"] = c["grid"]["ny"] = c["grid"]["nz"] = 5
    p = CA.parameters(c)
    states = ones(Int,5,5,5)
    octa = zeros(4,5,5,5)
    octa[4,3,3,3] = 0.0123
    temp = fill(1580.0,5,5,5)
    frozen = copy(octa)
    CA.updateDendriteLengths(states,octa,temp,p,0.001)
    @test octa == frozen # No mushy cells: retain existing size without growth.
    octa .= 0
    states[3,3,3] = 2
    CA.updateDendriteLengths(states,octa,temp,p,0.001)
    expected_increment = CA.dendriteTipVelocity(p.Tliquidus-1580,p)*0.001
    expected = zeros(5,5,5)
    expected[2:4,2:4,2:4] .= expected_increment
    expected[3,3,3] = 0
    @test octa[4,:,:,:] == expected
    @test count(>(0),octa[4,:,:,:]) == 26 # Includes diagonal edge/corner neighbours.
    @test octa[4,1,3,3] == 0 # Radius-two capture candidates cannot trigger growth.
    states .= 1
    states[1,1,1] = 2
    octa .= 0
    CA.updateDendriteLengths(states,octa,temp,p,0.001)
    @test count(>(0),octa[4,:,:,:]) == 7 # Clipped corner; no periodic wrapping.
    @test octa[4,5,1,1] == 0
end

if get(ENV,"AM_TEST_CCX","0") == "1"
    @testset "Real stock CalculiX thermal → CA integration" begin
        mktempdir() do dir
            result = run_calculix_demo(;resolution="quick",output_dir=dir,case="conduction")
            @test result.summary["solver_executed"]
            @test result.summary["solver_snapshots"] >= 2
            @test result.state.step == 100
            @test last(result.series.times_s) ≈ 0.002
            @test all(isfinite,result.state.previous_temperature)
            @test result.summary["finite_octahedra"]
            @test result.summary["octahedral_captures"] > 0
            # Bottom boundary remains fixed; a free interior node cools.
            @test result.series.temperature_K[1,end] ≈ 1500 atol=0.01
            @test result.series.temperature_K[23,end] < result.series.temperature_K[23,1]
            @test filesize(joinpath(dir,"final.vtk")) > 0
            replay = run_thermal_replay(;resolution="quick",dataset_dir=result.dataset_dir,output_dir=joinpath(dir,"replay"))
            @test !replay.summary["solver_executed"]
            @test result.state.ids == replay.state.ids
            @test result.state.octa == replay.state.octa
            @test result.state.grains == replay.state.grains
            @test result.rows == replay.rows
        end
    end
    @testset "Real layer heating, remelting and saved replay agreement" begin
        mktempdir() do dir
            result = run_calculix_demo(;resolution="quick",output_dir=dir,case="layer-heating")
            @test result.summary["solver_executed"]
            @test result.summary["ca_initialization"] == "solid-substrate"
            @test result.summary["melted_cells"] > 0
            @test result.summary["octahedral_captures"] > 0
            @test minimum(r.assigned_fraction for r in result.rows) < 0.5
            @test last(result.rows).assigned_fraction == 1.0
            replay = run_thermal_replay(;resolution="quick",dataset_dir=result.dataset_dir,output_dir=joinpath(dir,"replay"))
            @test result.state.ids == replay.state.ids
            @test result.state.octa == replay.state.octa
            @test result.state.grains == replay.state.grains
            @test result.rows == replay.rows
        end
    end
else
    @info "Real CalculiX test not requested; enable with AM_TEST_CCX=1 and CCX=/path/to/ccx"
end

@testset "Allocation guard before constructing CA arrays" begin
    @test_throws ArgumentError CA.check_grid(1000,1000,1000,8)
    @test_throws ArgumentError CA.check_grid(typemax(Int),typemax(Int),2,8)
    @test_throws ArgumentError CA.check_grid(10,10,10,8;max_cells=999)
    @test isnothing(CA.check_grid(10,10,10,8;max_cells=1000))
    c = CA.load_demo_config(CONFIG)
    c["grid"]["nx"] = 1_000_000
    @test_throws ArgumentError CA.parameters(c)
    mktempdir() do dir
        path = joinpath(dir,"oversized.toml")
        open(io->TOML.print(io,c),path,"w")
        @test_throws ArgumentError CA.load_demo_config(path)
    end
end
