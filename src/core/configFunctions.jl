# Description: Functions for handling reading, adjusting and creating new configuration files
using TOML
include("userInteractions.jl")

# Parameter struct declaration
if !isdefined(@__MODULE__, :Params)
Base.@kwdef struct Params
    Tliquidus::Int
    Tsolidus::Int
    tau::Float64
    
    lengthx::Float64
    lengthy::Float64
    lengthz::Float64

    nx::Int
    ny::Int
    nz::Int
    resolution::Int
    dx::Float64
    dt::Float64
    numGrains::Int
    minGrainSize::Float64

    D_L::Float64
    Gamma::Float64
    m_L::Float64
    k_0::Float64
    C_0::Float64

    dTsigma::Float64
    dTc::Float64
    N0::Float64

    configName::String
    filename::String
    multilayer::Bool
end

end # Allow configuration helpers to be included independently.

# Returns a Params struct with the configuration parameters. If the configuration file does not exist, it is created with default values
function getConfig(configFileName, lengthX, lengthY, lengthZ)
    endswith(configFileName, ".toml") ? nothing : configFileName = configFileName*".toml"
    # Check if configuration exist
    print("Searching for config file: ", configFileName)
    if (isfile(configFileName)) && (configFileHasAllFields(configFileName))
        p = from_toml(Params, configFileName)
        if (p.lengthx == lengthX && p.lengthy == lengthY && p.lengthz == lengthZ)
            refreshConfig(configFileName)
            return from_toml(Params, configFileName)
        end
    end
    createConfig(configFileName, lengthX, lengthY, lengthZ)
    wait_for_key("No existing config file with same dimensions was found. A config file with default parameters is created, please adjust it to your needs before pressing any button to continue with the calculations.")
    return refreshConfig(configFileName)
end

# Checks if the configuration file has all the fields required for the Params struct
function configFileHasAllFields(configFileName)
    open(configFileName) do file
        data = TOML.parse(file)
        for i in collect(fieldnames(Params))
            haskey(data, string(i)) ? nothing : return false # Checks all keys from toml to see if key i from Params is present, if not, return false
        end
        return true
    end 
end

# Creates a config file with default values, setting all dependent parameters to 0, and filling them with getDependentParams
function createConfig(configFileName, lengthX, lengthY, lengthZ)
    endswith(configFileName, ".toml") ? nothing : configFileName = configFileName*".toml"
    # Dezfoli, Lo & Raza (2021) for In718 parameters
    Tliquidus::Int = 1609
    Tsolidus::Int = 1533
    tau::Float64 = 1.0

    numGrains::Int = 100 #Not used anymore
    minGrainSize::Float64 = 0.01 #In mm, should be 0.002 for 2 μm

    D_L::Float64 = 1.6e-9
    Gamma::Float64 = 3.65e-7
    m_L::Float64 = -10.5
    k_0::Float64 = 0.48
    C_0::Float64 = 0.5

    dTsigma::Float64 = 6 
    dTc::Float64 = 12
    N0::Float64 = 1.2e6

    v = dendriteTipVelocity(Tliquidus-Tsolidus, D_L, Gamma, m_L, k_0, C_0)
    nx,ny,nz,resolution,dx,dt,filename = getDependentParams(lengthX, lengthY, lengthZ, minGrainSize, Tliquidus, Tsolidus, tau, v)
    filename = "crystalMeshRes"*string(resolution)
    p = Params(Tliquidus,Tsolidus,tau,
                  lengthX,lengthY,lengthZ, 
                  nx,ny,nz,resolution,dx,dt,numGrains,minGrainSize, 
                  D_L, Gamma, m_L, k_0, C_0, 
                  dTsigma, dTc, N0, 
                  configFileName, filename, 0)
    to_toml(configFileName,p)
    return p
end

# Refreshes the config file with loading new independendt parameters, and setting corresponding dependent parameters
function refreshConfig(configFileName)
    # Load config
    p = from_toml(Params, configFileName)
    v = dendriteTipVelocity(p.Tliquidus-p.Tsolidus, p.D_L, p.Gamma, p.m_L, p.k_0, p.C_0)
    # Get dependent parameters
    nx,ny,nz,resolution,dx,dt,filename = getDependentParams(p.lengthx, p.lengthy, p.lengthz, p.minGrainSize, p.Tliquidus, p.Tsolidus, p.tau, v)
    # Define new config, creating/replacing the dependent parameters with compatible ones
    newP = Params(p.Tliquidus,p.Tsolidus,p.tau,
                  p.lengthx,p.lengthy,p.lengthz, 
                  nx,ny,nz,resolution,dx,dt,p.numGrains,p.minGrainSize, 
                  p.D_L, p.Gamma, p.m_L, p.k_0, p.C_0, 
                  p.dTsigma, p.dTc, p.N0, 
                  configFileName, filename, p.multilayer)
    to_toml(configFileName,newP)
    return newP
end

# Sets dependent parameters 
function getDependentParams(lengthX, lengthY, lengthZ, minGrainSize, Tliquidus, Tsolidus, tau, v)
    nx,ny,nz = floor(Int, lengthX/minGrainSize), floor(Int, lengthY/minGrainSize), floor(Int, lengthZ/minGrainSize)
    resolution = max(nx,ny,nz)
    dx = min(lengthX/nx, min(lengthY/ny, lengthZ/nz))
    # dt = tau*dx/v #Eq 5.23 in (Beuerlein, 2023)
    dt = 3.541667e-05 # TODO: Replace this with a function that takes the steptime from ccx and combined with line above
    filename = "crystalMeshRes"*string(resolution)
    return nx,ny,nz,resolution,dx,dt, filename
end

# Returns the dendrite tip velocity, should be same as in crystalGrowth.jl
function dendriteTipVelocity(dT, D_L, Gamma, m_L, k_0, C_0) # TODO: prevent repeating this function, almost identical to the one in crystalGrowth.jl
    return (D_L / (5.51 * (pi^2) * ((-m_L * (1 - k_0))^1.5) * Gamma)) * ((dT.^2.5) / (C_0^1.5))
end


# Standalone configuration helpers use Julia's standard TOML library.
function from_toml(::Type{Params}, path::AbstractString)
    values = TOML.parsefile(path)
    return Params(; (Symbol(k) => v for (k, v) in values)...)
end

to_dict(p::Params) = Dict(string(k) => getfield(p, k) for k in fieldnames(Params))
function to_toml(path::AbstractString, p::Params)
    open(path, "w") do io
        TOML.print(io, to_dict(p); sorted=true)
    end
    return path
end
