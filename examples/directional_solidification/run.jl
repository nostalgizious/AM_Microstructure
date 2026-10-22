#!/usr/bin/env julia
using AMMicrostructure

function main(args)
    config = joinpath(@__DIR__, "config.toml")
    output = joinpath(@__DIR__, "..", "..", "output", "demo")
    i = 1
    while i <= length(args)
        if args[i] in ("--help", "-h")
            println("Usage: julia --project=. --threads=1 examples/directional_solidification/run.jl [--config FILE] [--output DIR]")
            return
        end
        args[i] in ("--config", "--output") || error("Unknown argument: $(args[i])")
        i < length(args) || error("Missing value for $(args[i])")
        args[i] == "--config" ? (config=args[i+1]) : (output=args[i+1])
        i += 2
    end
    AMMicrostructure.run_demo(;config_path=config,output_dir=output)
end

main(ARGS)
