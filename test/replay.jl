@testset "Fine CA preset preserves archived thermal geometry and duration" begin
    saved = CA.load_thermal_dataset(joinpath(@__DIR__,"..","data","layer_heating_demo"))
    before = deepcopy(saved.config)
    quick = CA.thermal_replay_config(saved.config,"quick")
    fine = CA.thermal_replay_config(saved.config,"fine")
    pq,pf = CA.parameters(quick),CA.parameters(fine)
    @test (pf.nx,pf.ny,pf.nz) == (128,24,80)
    @test pf.nx*pf.ny*pf.nz == 245760
    @test pf.dx == 0.0005
    @test pf.numGrains == 32
    @test fine["steps"] == 800
    @test pf.dt == 0.0025
    @test (pf.lengthx,pf.lengthy,pf.lengthz) == (pq.lengthx,pq.lengthy,pq.lengthz)
    @test fine["steps"]*pf.dt == quick["steps"]*pq.dt == 2.0
    @test fine["snapshot_step"]*pf.dt == quick["snapshot_step"]*pq.dt == 0.4
    @test saved.config == before
    @test quick["grid"] == saved.config["grid"]
    showcase = CA.thermal_replay_config(saved.config,"showcase")
    ps = CA.parameters(showcase)
    @test (ps.nx,ps.ny,ps.nz) == (384,72,240)
    @test ps.nx*ps.ny*ps.nz == 6635520
    @test ps.numGrains == 64
    @test ps.lengthx ≈ pq.lengthx
    @test showcase["max_cells"] == 7000000
    @test_throws ArgumentError CA.thermal_replay_config(saved.config,"unknown")
    # Mapping the finer cell centres still covers the same real FEM field.
    hot = map_temperatures(saved.mesh,saved.series,pf,0.4)
    @test size(hot) == (128,24,80)
    @test minimum(hot) >= 1500 && maximum(hot) <= 1800
end

@testset "Two real localized CalculiX heating passes and saved replay" begin
    dataset = joinpath(@__DIR__,"..","data","two_pass_heating_demo")
    saved = CA.load_thermal_dataset(dataset)
    @test saved.provenance["source_kind"] == "calculix"
    @test saved.provenance["thermal_case"] == "two-pass-heating"
    @test length(saved.mesh.node_ids) == 459
    @test length(saved.series.times_s) == 201
    @test last(saved.series.times_s) == 0.004
    top = saved.mesh.node_ids[:,2,end]
    hot1 = CA.nodal_temperatures(saved.series,0.0004)[top]
    hot2 = CA.nodal_temperatures(saved.series,0.0024)[top]
    @test argmax(hot1) == 7 # 24 μm first pulse centre.
    @test argmax(hot2) == 11 # 40 μm second pulse centre.
    @test maximum(hot1) ≈ 1800 atol=0.01
    @test maximum(hot2) ≈ 1800 atol=0.01
    @test first(hot1) < 1533 && last(hot2) < 1533
    @test maximum(CA.nodal_temperatures(saved.series,0.002)) < 1533
    @test maximum(CA.nodal_temperatures(saved.series,0.004)) < 1533
    mktempdir() do dir
        result = run_thermal_replay(;dataset_dir=dataset,output_dir=dir,resolution="quick")
        @test result.summary["heating_intervals_ms"] == [[0.0,0.4],[2.0,2.4]]
        @test result.state.step == 200
        @test result.summary["melted_cells"] > 0
        @test result.summary["octahedral_captures"] > 0
        @test result.summary["finite_octahedra"]
        @test isfile(joinpath(dir,"first_cycle.vtk")) && isfile(joinpath(dir,"second_peak.vtk"))
        @test minimum(r.assigned_fraction for r in result.rows if 0 < r.time_ms < 2) < 1
        @test minimum(r.assigned_fraction for r in result.rows if r.time_ms > 2) < 1
    end
end

@testset "Committed actual CalculiX dataset, solver-free replay" begin
    dataset = joinpath(@__DIR__,"..","data","calculix_demo")
    saved = CA.load_thermal_dataset(dataset)
    @test saved.provenance["source_kind"] == "calculix"
    @test saved.provenance["solver_version"] == "2.23"
    @test length(saved.series.times_s) == 101 # Includes explicit initial state.
    @test saved.series.times_s[end] ≈ 0.002
    @test saved.series.temperature_K[1,end] ≈ 1500 atol=0.01
    @test saved.series.temperature_K[23,end] < saved.series.temperature_K[23,1]
    mktempdir() do dir
        a = run_thermal_replay(;resolution="quick",dataset_dir=dataset,output_dir=joinpath(dir,"first"))
        b = run_thermal_replay(;resolution="quick",dataset_dir=dataset,output_dir=joinpath(dir,"second"))
        @test !a.summary["solver_executed"]
        @test a.state.ids == b.state.ids
        @test a.state.octa == b.state.octa
        @test a.state.grains == b.state.grains
        @test a.rows == b.rows
        @test a.summary["cells"] == 576
        @test a.summary["octahedral_captures"] > 0
        @test a.summary["assigned_fraction"] == 1.0
        @test filesize(joinpath(dir,"first","preview.svg")) > 0
        sparse_config = deepcopy(saved.config)
        sparse_config["steps"],sparse_config["slice_stride"] = 5,2
        sparse_dir = joinpath(dir,"sparse")
        sparse = CA.evolve_thermal_history(saved.mesh,saved.series,sparse_config,sparse_dir)
        @test sparse.summary["slice_frames"] == 4
        @test sparse.summary["slice_stride"] == 2
        exported = readlines(joinpath(sparse_dir,"slice_history.csv"))[2:end]
        @test unique(parse(Int,first(split(row,','))) for row in exported) == [0,2,4,5]
        @test length(sparse.rows) == 6 # Full scalar history includes unsaved slice steps.
        tampered = joinpath(dir,"tampered")
        cp(dataset,tampered)
        open(io->println(io,"changed"),joinpath(tampered,"thermal.dat"),"a")
        @test_throws ErrorException CA.load_thermal_dataset(tampered)
    end
end

@testset "Actual layer heating dataset melts and regrows the original CA" begin
    saved = CA.load_thermal_dataset(joinpath(@__DIR__,"..","data","layer_heating_demo"))
    @test saved.provenance["thermal_case"] == "layer-heating"
    @test length(saved.mesh.node_ids) == 135
    @test all(==(1500.0),saved.series.temperature_K[:,1])
    top = vec(saved.mesh.node_ids[:,:,end])
    hot = CA.nodal_temperatures(saved.series,0.0004)
    cold = CA.nodal_temperatures(saved.series,0.002)
    @test all(hot[top] .≈ 1800.0)
    @test all(cold[top] .< saved.config["material"]["solidus_K"])
    @test all(cold[top] .< hot[top])
    mktempdir() do dir
        result = run_thermal_replay(;resolution="quick",dataset_dir=joinpath(@__DIR__,"..","data","layer_heating_demo"),output_dir=dir)
        @test !result.summary["solver_executed"]
        @test result.summary["cells"] == 3840
        @test result.summary["melted_cells"] > 0
        @test result.summary["octahedral_captures"] > 0
        @test minimum(r.assigned_fraction for r in result.rows) < 0.5
        @test first(result.rows).assigned_fraction == last(result.rows).assigned_fraction == 1.0
        @test result.summary["finite_octahedra"]
        @test result.summary["resolution"] == "quick"
        @test result.summary["slice_frames"] == result.state.step+1
        @test result.summary["final_slice_grains"] == length(unique(result.state.ids[:,cld(result.state.params.ny,2),:]))
        @test TOML.parsefile(joinpath(dir,"config.toml"))["grid"] == saved.config["grid"]
        # Exported animation frames must describe the real CA state and timing.
        lines = readlines(joinpath(dir,"slice_history.csv"))
        fields = split.(lines[2:end],',')
        p = result.state.params
        plane_size = p.nx*p.nz
        @test length(fields) == (result.state.step+1)*plane_size
        @test sort(unique(parse(Int,r[1]) for r in fields)) == collect(0:result.state.step)
        @test all(parse(Int,r[4]) == cld(p.ny,2) for r in fields)
        final = fields[end-plane_size+1:end]
        @test [parse(Int,r[6]) for r in final] == vec(result.state.ids[:,cld(p.ny,2),:])
        @test [parse(Float64,r[8]) for r in final] ≈ vec(result.state.previous_temperature[:,cld(p.ny,2),:])
        @test parse(Float64,final[end][2]) ≈ result.state.time_ms
    end
end
