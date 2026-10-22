using AMMicrostructure: Params, getConfig, from_toml, to_dict, load_array, store_array, searchdir, wait_for_key
include("preinitializeCrystalMesh.jl")
include("findBaryCoords.jl")

using HDF5
using Printf

function findDims(meshVarsPtr, nodesPtr, numNodes) # C facing function, with nx,ny,nz returned
    nodes = zeros(Float64, 3*numNodes)
    load_array(nodesPtr, 3*numNodes, nodes)
    nodes = reshape(nodes, (3, :))
    lengthX, lengthY, lengthZ = maximum(nodes, dims=2) .- minimum(nodes, dims=2)
    configFileName = "crystalSimulationConfig.toml"
    p = getConfig(configFileName, lengthX, lengthY, lengthZ)
    nx, ny, nz, resolution, numGrains = p.nx, p.ny, p.nz, p.resolution, p.numGrains
    typeof(nx) != Int ? nx = convert(Int, nx) : nothing
    typeof(ny) != Int ? ny = convert(Int, ny) : nothing
    typeof(nz) != Int ? nz = convert(Int, nz) : nothing
    typeof(resolution) != Int ? resolution = convert(Int, resolution) : nothing
    typeof(numGrains) != Int ? numGrains = convert(Int, numGrains) : nothing
    prevHeight::Int = 0
    currHeight::Int = 0
    print(typeof(nx), " ", typeof(ny), " ", typeof(nz), " ", typeof(resolution), " ", typeof(numGrains), "\n")
    numGrains = floor(Int, nx*ny)   
    store_array(meshVarsPtr, 7,[nx, ny, nz, resolution, numGrains, prevHeight, currHeight])
end

function loadCrystalMesh(baryCoordsArrayPtr, elemArrayPtr, nodesPtr, numNodes, elemsPtr, numElems, ipkonPtr, ipkonLength, cellOctaPtr, cellIDArrPtr, grainArrPtr, nx, ny, nz, resolution, numGrains)
    dir = "/net/nas/results/240514_sf_CrystalMeshInterpolation/multiLayer/"
    println("Entered loadCrystalMesh")
    p = from_toml(Params, "crystalSimulationConfig.toml") #Has to match configFileName defined in findDims
    filename = p.filename*".h5"

    if p.nx != nx || p.ny != ny || p.nz != nz || p.resolution != resolution
        print("Parameters from the config file: ", p.nx, " ", p.ny, " ", p.nz, " ", p.resolution, "\n")
        print("Parameters from the mesh: ", nx, " ", ny, " ", nz, " ", resolution, "\n")
        error("The dimensions of the mesh does not match the dimensions of the config file. Please do not change the dimensions of the mesh during runtime.")
    end

    changed = 0
    if searchdir(dir, p.filename*"_grains.h5")  != []
        h5open(joinpath(dir, p.filename*"_grains.h5"), "r") do file
            cellIDArr = read(file, "cellIDArr")
            grainArr = read(file, "grainArr")
        end
    else 
        wait_for_key("No existing grain initialization for the crystal calculations was found. If this is not a prerun, restart. A new crystalMesh file will now be created, based on the supplied config file.")
        changed = 1
    end
    if searchdir(dir, p.filename*"_geom.h5")  != []
        h5open(joinpath(dir,p.filename*"_geom.h5"), "r") do file
            baryCoordsArray = read(file, "baryCoordsArray")
            elemArray = read(file, "elemArray")
        end
    else
        wait_for_key("No existing mesh interpolation for the crystal calculations was found. If this is not a prerun, restart. A new crystalMesh file will now be created, based on the supplied config file.")
        baryCoordsArray, elemArray = calcBaryCoords(baryCoordsArrayPtr, elemArrayPtr, nodesPtr, numNodes, elemsPtr, numElems, ipkonPtr, ipkonLength, p) #Storing to C is done within the function
        changed = 1
    end
    if(changed == 1)
    end
    
    grainConfig = h5read(joinpath(dir,p.filename*"_grains.h5"), "config")
    geomConfig = h5read(joinpath(dir,p.filename*"_geom.h5"), "config")
    paramFileConfig = @sprintf("%s",(to_dict(p)))
    if(grainConfig != paramFileConfig || geomConfig != paramFileConfig)
        wait_for_key("The existing initial mesh does not exactly match the config file. Ensure that no essential changes have been made.")
        print("Parameters in the saved initial geometric mesh: \n", geomConfig, "\n")
        print("Parameters in the saved initial grains: \n", grainConfig, "\n")
        print("Parameters in the config file: \n", paramFileConfig, "\n")
        wait_for_key("Please double check.")
    end

    cellOctaArr = zeros(Float64, 4*nx*ny*nz)
    load_array(cellOctaPtr, 4*nx*ny*nz, cellOctaArr)
    cellOctaArr = reshape(cellOctaArr, (4, nx, ny, nz))
    Threads.@threads for i = 1:nx for j = 1:ny for k = 1:nz
                cellOctaArr[1:3,i,j,k] = [(i)/(p.nx+1)*p.lengthx, (j)/(p.ny+1)*p.lengthy, (k)/(p.nz+1)*p.lengthz]
            end end end
    store_array(cellOctaPtr, cellOctaArr)
    store_array(cellIDArrPtr, cellIDArr)
    store_array(grainArrPtr, grainArr)
    store_array(baryCoordsArrayPtr, baryCoordsArray)
    store_array(elemArrayPtr, elemArray)
end

function correctCrystalMeshExists(p, filename)
    if searchdir(".", filename)  != []
        h5open(filename, "r") do file
            stringedConfig = @sprintf("%s",(to_dict(p)))
            if haskey(file, "config") && read(file,"config")==stringedConfig && haskey(file, "cellStateArr") && haskey(file, "cellIDArr") && haskey(file, "grainArr")
                return true
            else
                return false
            end
        end
    else
        return false
    end
end

function correctMeshInterpolationExists(p, filename)
    if searchdir(".", filename)  != []
        h5open(filename, "r") do file
            stringedConfig = @sprintf("%s",(to_dict(p)))
            if haskey(file, "config") && read(file,"config")==stringedConfig && haskey(file, "baryCoordsArray") && haskey(file, "elemArray")
                return true
            end
            return false
        end
    else
        return false
    end
end
