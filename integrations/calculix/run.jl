#!/usr/bin/env julia
using AMMicrostructure

function main(args)
    ccx = get(ENV,"CCX","ccx")
    output = joinpath(@__DIR__,"..","..","output","calculix")
    case = "two-pass-heating"
    resolution = "showcase"
    i=1
    while i <= length(args)
        if args[i] in ("--help","-h")
            println("Usage: julia --project=. --threads=8,0 integrations/calculix/run.jl [--ccx EXECUTABLE] [--output DIR] [--case two-pass-heating|layer-heating|conduction] [--resolution showcase|fine|quick]")
            return
        end
        args[i] in ("--ccx","--output","--case","--resolution") || error("Unknown argument: $(args[i])")
        i < length(args) || error("Missing argument value")
        args[i] == "--ccx" ? (ccx=args[i+1]) : args[i] == "--case" ? (case=args[i+1]) : args[i] == "--resolution" ? (resolution=args[i+1]) : (output=args[i+1])
        i+=2
    end
    run_calculix_demo(;ccx,output_dir=output,case,resolution)
end
main(ARGS)
