using Test, LinearAlgebra, Random, TOML
using AMMicrostructure
const CA = AMMicrostructure
const CONFIG = joinpath(@__DIR__, "..", "examples", "directional_solidification", "config.toml")

@testset "Euler rotation and neighborhoods" begin
    v = [0.2,-0.7,1.3]
    angles = [23.0,-34.0,17.0]
    rotated = CA.rotVector(v, angles...)
    @test norm(rotated) ≈ norm(v)
    @test CA.rotVector(rotated, -reverse(angles)...) ≈ v atol=1e-14
    neighbors = CA.getMooreN(2)
    @test length(unique(neighbors)) == 124
    @test CartesianIndex(0,0,0) ∉ neighbors
    @test all(maximum(abs,Tuple(n)) <= 2 for n in neighbors)
end

@testset "Historical kinetics and Gaussian quadrature" begin
    p = CA.parameters(CA.load_demo_config(CONFIG))
    @test CA.dendriteTipVelocity(0.0,p) == 0
    @test CA.dendriteTipVelocity(20.0,p) ≈ 0.015983680761234173
    @test CA.dendriteTipVelocity(20.0,p) > CA.dendriteTipVelocity(10.0,p) > 0
    @test CA.nucleationIntegral(2.0,3.0,p) == 0
    # Full Gaussian mass includes sigma: retain original law's normalization.
    @test CA.nucleationIntegral(p.dTc+12p.dTsigma,p.dTc-12p.dTsigma,p) ≈ p.N0*p.dTsigma*p.dx^3 rtol=1e-8
    @test CA.nucleationIntegral(p.dTc,p.dTc-12p.dTsigma,p) ≈ 0.5p.N0*p.dTsigma*p.dx^3 rtol=1e-8
    a,b,c = 2.0,11.0,24.0
    @test CA.nucleationIntegral(c,a,p) ≈ CA.nucleationIntegral(b,a,p)+CA.nucleationIntegral(c,b,p) rtol=1e-6
end

@testset "Physical centers, melting and capture geometry" begin
    c = CA.load_demo_config(CONFIG)
    c["grid"]["nx"] = c["grid"]["ny"] = c["grid"]["nz"] = 3
    p = CA.parameters(c)
    ids, states, octa = zeros(Int,3,3,3), fill(2,3,3,3), zeros(4,3,3,3)
    grains = zeros(3,27)
    source = CartesianIndex(2,2,2)
    states[source], ids[source] = 1, 1
    @test CA.cellCoords(1,1,1,p,octa) == fill(p.dx/2,3)
    @test CA.cellCoords(3,3,3,p,octa) ≈ fill(2.5p.dx,3)
    grow!() = CA.updateCellGrowth(states,ids,octa,grains,p,CA.getMooreN(2),CA.toCartesianInd(CA.N1Ind),CA.OctoNormVec,3)
    @test first(grow!()) == 0
    # Previously the opposite half-space alone falsely captured distant cells.
    octa[4,source] = 0.99p.dx/sqrt(3)
    @test first(grow!()) == 0
    octa[4,source] = 1.01p.dx/sqrt(3)
    captures, _, faces = grow!()
    @test captures == sum(faces) == 6
    @test count(==(1),states) == 7 # Captured cells cannot capture again this sweep.
    @test ids[1,1,1] == 0
    @test all(isfinite,octa)
    @test all(octa[4,:,:,:] .>= 0)
    temp = fill(Float64(p.Tliquidus),size(ids))
    CA.updateCellStates(states,ids,octa,temp,temp,p,Dict(),3)
    @test all(==(3),states) && all(iszero,ids) && all(iszero,octa)
    cold = fill(Float64(p.Tsolidus-10),size(ids))
    CA.updateCellStates(states,ids,octa,temp,cold,p,Dict(),3)
    @test all(==(2),states) # Fast cooling must not strand unassigned liquid.
end

@testset "Decentered child rotation covariance" begin
    p = CA.parameters(CA.load_demo_config(CONFIG))
    origin = [0.01,0.02,0.03]
    offset = [0.002,0.0004,0.0002]
    normal = fill(1/sqrt(3),3)
    L = 0.002
    c0,l0 = CA.calcCentreOfEngulfedCell(zeros(3),origin,origin+offset,L,normal,normal,p)
    angles = [23.0,-34.0,17.0]
    nrot = CA.rotVector(normal,angles...)
    c1,l1 = CA.calcCentreOfEngulfedCell(angles,origin,origin+CA.rotVector(offset,angles...),L,nrot,normal,p)
    @test l1 ≈ l0 atol=1e-14
    @test c1-origin ≈ CA.rotVector(c0-origin,angles...) atol=1e-14
    @test l1 > 0 && all(isfinite,c1)
end

@testset "Standalone run, repeatability and output" begin
    mktempdir() do dir
        c = CA.load_demo_config(CONFIG)
        c["grid"]["nx"], c["grid"]["ny"], c["grid"]["nz"] = 10,3,8
        c["grid"]["substrate_grains"] = 4
        c["steps"] = 60
        cfg = joinpath(dir,"config.toml")
        open(io->TOML.print(io,c),cfg,"w")
        a = CA.run_demo(;config_path=cfg,output_dir=joinpath(dir,"a"))
        b = CA.run_demo(;config_path=cfg,output_dir=joinpath(dir,"b"))
        @test a.ids == b.ids
        @test a.octa == b.octa
        @test a.rows == b.rows
        @test a.summary["octahedral_captures"] > 0
        @test a.summary["assigned_fraction"] > first(a.rows).fraction
        @test a.summary["finite_octahedra"] && a.summary["valid_grain_ids"]
        @test all(diff([r.fraction for r in a.rows]) .>= 0)
        @test count(>(0),a.ids) == 30 + a.summary["nucleated_grains"] + a.summary["octahedral_captures"]
        for name in ("initial.vtk","midpoint.vtk","final.vtk","preview.svg","history.csv","grains.csv","summary.toml","config.toml")
            @test filesize(joinpath(dir,"a",name)) > 0
        end
        @test occursin("CELL_DATA 240",read(joinpath(dir,"a","final.vtk"),String))
        @test length(readlines(joinpath(dir,"a","history.csv"))) == 62
        p = CA.parameters(c)
        path = joinpath(dir,"params.toml")
        CA.to_toml(path,p)
        @test CA.to_dict(CA.from_toml(CA.Params,path)) == CA.to_dict(p)
        @test CA.configFileHasAllFields(path)
        malformed = CA.to_dict(p)
        malformed["not_nx"] = pop!(malformed,"nx")
        open(io->TOML.print(io,malformed),path,"w")
        @test !CA.configFileHasAllFields(path)
        retained = CA.to_dict(p)
        retained["multilayer"] = true
        open(io->TOML.print(io,retained),path,"w")
        @test CA.refreshConfig(path).multilayer
        created = CA.createConfig(joinpath(dir,"new.toml"),0.1,0.1,0.1)
        @test created.nx == created.ny == created.nz == 10
    end
end

include("fem.jl")

include("replay.jl")
