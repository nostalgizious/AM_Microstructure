#!/usr/bin/env julia
using AMMicrostructure

function main(args)
    dataset = joinpath(@__DIR__,"..","..","data","two_pass_heating_demo")
    output = joinpath(@__DIR__,"..","..","output","fem_replay")
    resolution = "showcase"
    i=1
    while i <= length(args)
        if args[i] in ("--help","-h")
            println("Usage: julia --project=. --threads=8,0 examples/fem_replay/run.jl [--dataset DIR] [--output DIR] [--resolution showcase|fine|quick]")
            return
        end
        args[i] in ("--dataset","--output","--resolution") || error("Unknown argument: $(args[i])")
        i < length(args) || error("Missing argument value")
        args[i] == "--dataset" ? (dataset=args[i+1]) : args[i] == "--resolution" ? (resolution=args[i+1]) : (output=args[i+1])
        i+=2
    end
    run_thermal_replay(;dataset_dir=dataset,output_dir=output,resolution)
end
main(ARGS)
