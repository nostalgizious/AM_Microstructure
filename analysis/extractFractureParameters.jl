using Statistics  # For mean, std, minimum, and maximum functions
using StatsBase
using Random
using HDF5
using Plots
using LinearAlgebra
using LaTeXStrings
using LsqFit

using AMMicrostructure: searchdir

macro Name(arg)
    string(arg)
 end

# resdir = "/net/nas/results/"
# projname = "DebuggingCrystalGrowth"
# simname = "Simulation"
# filename = ""
# resname = getResName(resdir, projname, simname, filename)
# # resname = "/net/nas/results/240430_sf_DebuggingCrystalGrowth/240430_1636_sf_Simulation"
# p = from_toml(Params, joinpath(resname, "config.toml"))

# nx, ny, nz = p.nx, p.ny, p.nz

# filename = "Results.h5"
# # resname = joinpath(resname, filename)
# h5name = joinpath(resname, filename)
# cellIDArr = h5read(h5name, "cellIDArr")
# grainArr = h5read(h5name, "grainArr")

# # trim of edges of mesh
# cellIDArr = cellIDArr[floor(Int, 0.1*p.nx):floor(Int, 0.9*p.nx), floor(Int, 0.1*p.ny):floor(Int, 0.9*p.ny), 1:p.nz-1]

# # grainArr = grainArr[:, 1:count(grainArr[1, :] .!= 0)]  # Remove zeros from the array
# numGrains = count(grainArr[1,:] .!= 0)

function trimCellID!(cellIDArr)
    height = 0
    nx, ny, nz = size(cellIDArr)
    for i=1:nz
        count(cellIDArr[:,:,i] .==0)/(nx*ny) < 0.95 ? (height += 1) : break
    end
    # cellIDArr = cellIDArr[:,:,1:height]
    return height
end

# Sets all preinitialized areas to zero, and returns a cropped area containing only a bounding box including the recrystallized cells
function sLtrimCellID(cellIDArr, cellStateArr)
    cellIDArr[cellStateArr .== 0] .= 0
    nx, ny, nz = size(cellIDArr)
    minX, minY, minZ = nx, ny, nz
    maxX, maxY, maxZ = 1, 1, 1
    count(cellIDArr .!= 0) == 0 ? (return Array{eltype(cellIDArr)}(undef, 0, 0, 0)) : nothing
    R = CartesianIndices(cellIDArr)
    Threads.@threads for ind in R
        cellIDArr[ind] .== 0 ? continue : nothing
        minX = min(minX, ind[1])
        minY = min(minY, ind[2])
        minZ = min(minZ, ind[3])
        maxX = max(maxX, ind[1])
        maxY = max(maxY, ind[2])
        maxZ = max(maxZ, ind[3])
    end
    # minX, maxX, minY, maxY, minZ, maxZ
    return cellIDArr[minX:maxX, minY:maxY, minZ:maxZ]
end

function trimGrainArr!(grainArr, cellIDArr)
    grainArr = grainArr[:, 1:maximum(cellIDArr)]
end


# Function to calculate the aspect ratio of each grain
function getAreaVolume(cellIDArr)
    cellIDArr = cellIDArr[cellIDArr .!= 0]
    grains = unique(cellIDArr)
    volumes = zeros(2,length(grains))
    Threads.@threads for gnum in 1:length(grains)
        g = grains[gnum]
        # Find indices where the current grain ID is located
        # indices = findall(x -> x == g, cellIDArr)
        indices = findall(cellIDArr .== g)
        volumes[1,gnum] = g
        volumes[2,gnum] = length(indices)
    end
    return volumes.*(p.dx*1e3)^3
end

# function getEllipse(cellIDSlice, dx)
#     grains = unique(cellIDArr)
#     ellipsoidAxisLengths = zeros(length(grains), 2)
#     ellipsoidAxisDirs = zeros(length(grains), 2, 2)
#     Threads.@threads for g in 1:length(grains)
#         # Find indices where the current grain ID is located
#         # indices = findall(x -> x == g, cellIDArr)
#         indices = findall(cellIDArr .== grains[g])
#         length(indices) == 0 ? (ellipsoidAxisLengths[g,1] =0; ellipsoidAxisDirs .= 0;continue) : n = length(indices) # Set grains without associated cells to -1
#         n == 1 ? (ellipsoidAxisLengths[g,1] = -2; ellipsoidAxisDirs .= 0; continue) : nothing # Set grains containing only a single cell to -2
#         centre = [sum(indices)[i] for i=1:2]./n # in units p.dx, not m or μm
#         gyrTensor = zeros(Float64, 2,2)
#         [(r = [ind[1], ind[2]] .-centre; gyrTensor += r*r') for ind in indices]./n
#         F = eigen(gyrTensor)
#         ellipsoidAxisLengths[g, :] = F.values
#         ellipsoidAxisDirs[g, :, :] = F.vectors
#     end
#     return ellipsoidAxisLengths, ellipsoidAxisDirs
# end

function getEllipseSimple(cellIDArr, dx)
    grains = unique(cellIDArr)
    ellipsoidAxisLengths = zeros(Float64, length(grains), 2)
    ellipsoidAxisDirs = zeros(Float64, length(grains), 2, 2)
    areas = zeros(Float64, length(grains))
    counters = zeros(Int32, Threads.nthreads())
    Threads.@threads for g in 1:length(grains)
        counters[Threads.threadid()] += 1
        # Find indices where the current grain ID is located
        # indices = findall(x -> x == g, cellIDArr)
        indices = findall(cellIDArr .== grains[g])
        length(indices) == 0 ? (ellipsoidAxisLengths[g,1] = 0; continue) : n = length(indices) # Set grains without associated cells to -1
        n == 1 ? (ellipsoidAxisLengths[g,1] = 0; continue) : nothing # Set grains containing only a single cell to -2
        
        centre = [sum(indices)[i] for i=1:2]./n # in units p.dx, not m or μm
        gyrTensor = zeros(Float64, 2,2)
        [(r = [ind[1], ind[2]] .-centre; gyrTensor += r*r') for ind in indices]./n
        F = eigen(gyrTensor)
        ellipsoidAxisLengths[g, :] = F.values #* (dx*1e3)
        ellipsoidAxisDirs[g, :, :] = F.vectors

        areas[g] = length(indices)
        # counters[Threads.threadid()]%1000 == (Threads.threadid()*100) ? println("Percent completed: ", sum(counters)/length(grains)*100) : nothing
    end
    areas = areas.*(dx*1e3)^2 #Gives result in μm³ if dx is specified in mm
    # ellipsoidAxisLengths = ellipsoidAxisLengths*(dx*1e3) # Gives result in μm if dx is specified in mm
    return ellipsoidAxisLengths, ellipsoidAxisDirs, areas
end

function getEllipse(cellIDArr, dx)
    grains = unique(cellIDArr)
    ellipsoidAxisLengths = zeros(Float64, length(grains), 2)
    ellipsoidAxisDirs = zeros(Float64, length(grains), 2, 2)
    areas = zeros(Float64, length(grains))
    counters = zeros(Int32, Threads.nthreads())
    Threads.@threads for g in 1:length(grains)
        counters[Threads.threadid()] += 1
        # Find indices where the current grain ID is located
        # indices = findall(x -> x == g, cellIDArr)
        indices = findall(cellIDArr .== grains[g])
        length(indices) == 0 ? (ellipsoidAxisLengths[g,1] = 0; continue) : n = length(indices) # Set grains without associated cells to -1
        n == 1 ? (ellipsoidAxisLengths[g,1] = 0; continue) : nothing # Set grains containing only a single cell to -2
        # centres = [sum(indices)[i] for i=1:2]./n # in units p.dx, not m or μm

        centres = zeros(n,2)
        for i=1:n #Attempt at removing outliers, drastically reduced performance
            centres[i,:] = [sum(indices[1:end .!= i])[j] for j=1:2]
        end
        centre = median.(eachcol(centres))./(n-1)
        rs = [([ind[1], ind[2]] .-centre) for ind in indices]
        norms = norm.(eachrow(rs))
        rs = rs[norms .< percentile(norms, 80), :]
        gyrTensor = zeros(Float64, 2,2)
        # println(rs)
        for r in rs
            # print("   ", r)
            gyrTensor += r*r'
        end
        gyrTensor = gyrTensor ./ (length(norms))

        # [(gyrTensor += r*r') for ind in indices]./n
        F = eigen(gyrTensor)
        ellipsoidAxisLengths[g, :] = F.values #* (dx*1e3)
        ellipsoidAxisDirs[g, :, :] = F.vectors
        areas[g] = length(indices)
        counters[Threads.threadid()]%1000 == (Threads.threadid()*100) ? println("Percent completed: ", sum(counters)/length(grains)*100) : nothing
    end
    areas = areas.*(dx*1e3)^2 #Gives result in μm³ if dx is specified in mm
    ellipsoidAxisLengths = ellipsoidAxisLengths*(dx*1e3) # Gives result in μm if dx is specified in mm
    return ellipsoidAxisLengths, ellipsoidAxisDirs, areas
end

function getEllipsoid(cellIDArr, dx)
    grains = unique(cellIDArr)
    ellipsoidAxisLengths = zeros(Float64, length(grains), 3)
    ellipsoidAxisDirs = zeros(Float64, length(grains), 3, 3)
    volumes = zeros(Float64, length(grains))
    counters = zeros(Int32, Threads.nthreads())
    Threads.@threads for g in 1:length(grains)
        counters[Threads.threadid()] += 1
        # Find indices where the current grain ID is located
        # indices = findall(x -> x == g, cellIDArr)
        indices = findall(cellIDArr .== grains[g])
        length(indices) == 0 ? (ellipsoidAxisLengths[g,1] = 0; continue) : n = length(indices) # Set grains without associated cells to -1
        n == 1 ? (ellipsoidAxisLengths[g,1] = 0; continue) : nothing # Set grains containing only a single cell to -2
        centre = [sum(indices)[i] for i=1:3]./n # in units p.dx, not m or μm
        gyrTensor = zeros(Float64, 3,3)
        [(r = [ind[1], ind[2], ind[3]] .-centre; gyrTensor += r*r') for ind in indices]./n
        F = eigen(gyrTensor)
        ellipsoidAxisLengths[g, :] = F.values #* (dx*1e3)
        ellipsoidAxisDirs[g, :, :] = F.vectors
        volumes[g] = length(indices)
        counters[Threads.threadid()]%1000 == (Threads.threadid()*100) ? println("Percent completed: ", sum(counters)/length(grains)*100) : nothing
    end
    volumes = volumes.*(dx*1e3)^3 #Gives result in μm³ if dx is specified in mm
    ellipsoidAxisLengths = ellipsoidAxisLengths*(dx*1e3) # Gives result in μm if dx is specified in mm
    return ellipsoidAxisLengths, ellipsoidAxisDirs, volumes
end

function getAspectRatios(eL,dx)
    # eL[eL[:,1] .> -1, :]
    # eL[eL .== 0] .= 0eL
    aspectRatio = zeros(length(eL[:,1]))
    Threads.@threads for i=1:length(eL[:,1])
        dmax = maximum(eL[i,:])
        dmin = max(minimum(eL[i,:]), dx)
        aspectRatio[i] = dmax/dmin
    end
    # aspectRatio = [maximum(eL[i,:]) / minimum(eL[i,:]) for i in 1:length(eL[:,1])]
    # aspectRatio[aspectRatio .== Inf .|| aspectRatio .== NaN] .= 0
    aspectRatio[aspectRatio .> 1e3 .|| aspectRatio .< 0] .= 0
    return aspectRatio
end

# vxmid = getAreaVolume(cellIDArr[floor(Int, p.nx/2),:,:], grainArr)
# vymid = getAreaVolume(cellIDArr[:,floor(Int, p.ny/2),:], grainArr)
# vzmid = getAreaVolume(cellIDArr[:,:,floor(Int, p.nz/2)], grainArr)

binrange = range(0,100,10)
# histogram(vxmid[vxmid .!= 0], bins=binrange, label="X", xlabel="Volume (um^3)", ylabel="Frequency", title="Volume distribution along X-axis")
# histogram(vymid[vymid .!= 0], bins=binrange, label="Y", reuse = false)
# histogram(vzmid[vzmid .!= 0], bins=binrange, label="Z", reuse = false)


# Calculate aspect ratios
# aspect_ratios = calculateAspectRatios(cellIDArr, grainArr)

# println("Aspect ratios of grains: $aspect_ratios")

function loadRelCells(h5n)
    cellIDArr = h5read(h5n,"cellIDArr")
    cellStateArr = h5read(h5n,"cellStateArr")
    relCells = cellIDArr[cellStateArr .== 1]
    return relCells
end

function loadRelGrains(h5n, cellIDArr)
    grainArr = h5read(h5n,"grainArr")
    return grainArr[:,unique(cellIDArr)]
end

function getSlicedStats(cellIDArr,dx, name)
    # cellIDArr = cellIDArr[:,:,1:trimCellID!(cellIDArr)]
    nx, ny, nz = size(cellIDArr)
    testx025 = floor(Int, nx/4)
    testx050 = floor(Int, nx/2)
    testy025 = floor(Int, ny/4)
    testy050 = floor(Int, ny/2)
    testz025 = floor(Int, 0.25*nz)
    testz050 = floor(Int, 0.50*nz)
    testz075 = floor(Int, 0.75*nz)
    slicedStats =Dict{String, Any}()
    slicedStats["ex025"] =      getEllipse(cellIDArr[testx025,:,:],dx)
    slicedStats["ex050"] =      getEllipse(cellIDArr[testx050,:,:],dx)
    slicedStats["ey025"] =      getEllipse(cellIDArr[:,testy025,:],dx)
    slicedStats["ey050"] =      getEllipse(cellIDArr[:,testy050,:],dx)
    slicedStats["ex025lower"] = getEllipse(cellIDArr[testx025,:,1:testz050],dx)
    slicedStats["ex050lower"] = getEllipse(cellIDArr[testx050,:,1:testz050],dx)
    slicedStats["ey025lower"] = getEllipse(cellIDArr[:,testy025,1:testz050],dx)
    slicedStats["ey050lower"] = getEllipse(cellIDArr[:,testy050,1:testz050],dx)
    slicedStats["ex025upper"] = getEllipse(cellIDArr[testx025,:,testz050:nz],dx)
    slicedStats["ex050upper"] = getEllipse(cellIDArr[testx050,:,testz050:nz],dx)
    slicedStats["ey025upper"] = getEllipse(cellIDArr[:,testy025,testz050:nz],dx)
    slicedStats["ey050upper"] = getEllipse(cellIDArr[:,testy050,testz050:nz],dx)
    slicedStats["ez025"] =      getEllipse(cellIDArr[:,:,testz025],dx)
    slicedStats["ez050"] =      getEllipse(cellIDArr[:,:,testz050],dx)
    slicedStats["ez075"] =      getEllipse(cellIDArr[:,:,testz075],dx)
    
    constKeys = deepcopy(keys(slicedStats))
    for key in constKeys
        localDict = Dict{String, Any}()
        i = slicedStats[key]
        eL, eDir, area = i[1], i[2], i[3]
        invChecker = [prod(eL[i,:]) for i in 1:length(eL[:,1])] # helper vector to determine indices of proper values
        indsTrimInv = findall(invChecker .> 0 .&& invChecker .!= NaN .&& invChecker .< 10^5) # find indices of valid entries
        indsTrimOutliers = findall(invChecker .> percentile(invChecker[indsTrimInv], 10) .&& invChecker .< percentile(invChecker[indsTrimInv], 99)) #find indices of middle 98%
        eLen = eL[indsTrimOutliers,:]
        eD = [norm(eLen[i,:]) for i in 1:length(eLen[:,1])]
        areaRaw = area *  (dx*1000)^2
        area= areaRaw[indsTrimOutliers]
        eArea = pi.*[prod(eLen[i,:]) for i in 1:length(eLen[:,1])]*(dx*1000)^2 #Scale to μm²
        AR = getAspectRatios(eLen, dx)
        # eVol[eVol .== NaN .|| eVol .== Inf .|| eVol .== -1] .= 0
        localDict["ellipseLengths"] = eLen
        localDict["grainRadius"] = eD
        localDict["ellipseArea"] = eArea
        localDict["aspectRatio"] = AR
        localDict["ellipseDirections"] = eDir
        localDict["ellipseLengthsRaw"] = eL
        localDict["countedArea"] = area
        localDict["ratioInvalid"] = 1 - length(indsTrimInv)/length(eL[:,1])
        localDict["ratioTrimmer"] = 1 - length(indsTrimOutliers)/length(indsTrimInv)
        localDict["areaWeightedAspectRatio"] = sum(AR .* area) / sum(area)
        localDict["elAreaWeightedAspectRatio"] = sum(AR .* eArea) / sum(eArea)
        localDict["sliceName"] = key
        localDict["rawArea"] = areaRaw
        # slicedStats[key] = nothing
        slicedStats[key] = localDict
        # slicedStats["derived_from_"*key] = localDict
    end
    # ax = getAreaVolume(cellIDArr[testx,:,:])
    # ay = getAreaVolume(cellIDArr[:,testy,:])
    # az1 = getAreaVolume(cellIDArr[:,:,testz1])
    # az2 = getAreaVolume(cellIDArr[:,:,testz2])
    # az3 = getAreaVolume(cellIDArr[:,:,testz3])
    slicedStats["simulatioName"] = name
    return slicedStats
end

function getSlicedGrainSize(cellIDArr)
    nx, ny, nz = size(cellIDArr)
    nz = trimCellID!(cellIDArr)
    startx, starty = floor(Int, 0.2*nx), floor(Int, 0.2*ny)
    stopx, stopy, stopz = ceil(Int, 0.8*nx), ceil(Int, 0.8*ny),  ceil(Int, 0.8*nz)
    redCellID = cellIDArr[startx:stopx, starty:stopy, 1:stopz]
    dist = zeros(3, stopz)
    Threads.@threads for i=1:stopz
        tmp = getAreaVolume(redCellID[:,:,i])
        areas = tmp[2,:]
        dist[1,i] = mean(areas)
        dist[2,i] = std(areas)
        dist[3,i] = length(areas)
    end
    return dist
end

function getVolumeWeighterDirections() #TODO: Transform ellipse dirs into local coordinate system of grain to show directional preference. Volume weighted as the ones oriented closely aligned with bd should be larger. idk, lets check

end

# pf = getResName("/net/nas/results/", "DebuggingCrystalGrowth", "Simulation", "config.toml")
# h5n = getResName("/net/nas/results/", "DebuggingCrystalGrowth", "Simulation", "Results.h5")
# p = from_toml(Params, pf)
# cellIDArr = h5read(h5n,"cellIDArr")

# as = getStats(cellIDArr)
# for i=1:lastindex(as)
#     h = histogram(as[i][2,:], bins=binrange)
#     display(h)
# end



function analyzeConv(outdir, namefragment="", arrayName="")
    files = searchdir(outdir,namefragment)
    length(files) > 1 ? error("Specify more of the name, several found") : nothing
    length(files) == 0 ? error("No matching files found") : nothing
    i = files[1]
    # simname = (replace(i, ".h5"=>""))*"_ellipses.h5"
    outpath = joinpath(outdir, i)
    # !haskey(outpath, arrayName) ? (println("Only these keys available:", h5open(outpath) do file return keys(file) end); error("Enter valid key.")) : nothing
    return h5open(outpath) do file 
        haskey(file, arrayName) ? (return read(file, arrayName)) :  (println("Only these keys available:", keys(file)); error("Enter valid key."))
    end
end


function compareGrainSizes(ml=false)
    if ml==true
        outdir = mloutdir
    else
        outdir = sloutdir #only for sl
    end
    meanGrainSizes = Vector{Float64}()
    meanGrainVolumes = Vector{Float64}()
    # binrange = range(0, 800, 10)
    for i in searchdir(outdir,"")
        println(i)
        d = analyzeConv(outdir, i, "countedVolume")
        # d = d[d .< percentile(d, 99) .&& d .> percentile(d, 1)]
        println(length(d))
        h = histogram(d, bins=binrange)
        display(h)
        push!(meanGrainSizes, mean(d))
    end
    for i in 1:length(searchdir(outdir,""))
        println(searchdir(outdir,"")[i], ": ", meanGrainSizes[i])
    end
    return meanGrainSizes
end


function checkConvCriteria(ml=false)
    if ml==true
        outdir = mloutdir
    else
        outdir = sloutdir #only for sl
    end
    vwars = Vector{Float64}()
    binrange = range(0.1, 20, 20)
    for i in searchdir(outdir,"")
        println(i)
        vwar = analyzeConv(outdir, i, "volumeWeightedAspectRatio")
        vwar = vwar[vwar .< percentile(vwar, 95) .&& vwar .> percentile(vwar, 5)]
        println(length(vwar))
        println(sum(vwar))
        h = histogram(vwar)
        display(h)
        push!(vwars, median(vwar)*length(vwar))
    end
    return vwars
end

function get3DSizes(ml=false)
    if ml==true
        outdir = mloutdir
    else
        outdir = sloutdir #only for sl
    end
    vwars = Vector{Float64}()
    binrange = range(0.1, 20, 20)
    for i in searchdir(outdir,"")
        println(i)
        vwar = analyzeConv(outdir, i, "countedVolume")
        println(length(vwar))
        vwar = vwar[vwar .> 0 .&& vwar .< 10^5]
        vwar = vwar[vwar .< percentile(vwar, 99) .&& vwar .> percentile(vwar, 10)]
        println(length(vwar))
        # println(sum(vwar))
        h = histogram(vwar)
        display(h)
        push!(vwars, mean(vwar))
    end
    return vwars
end

function runConvAnalysis(namefragment="", ml=false)
    if ml == true
        indir = mlindir
        outdir = mloutdir
    else
        indir = slindir
        outdir = sloutdir
    end
for i in searchdir(indir,namefragment)
    inpath = joinpath(indir, i)
    simname = (replace(i, ".h5"=>""))*"_ellipses.h5"
    outpath = joinpath(outdir, simname)
    println("Saving to :", outpath)
    isfile(outpath) ? continue : nothing
    h5write(outpath, "SimulationName", simname)
    cellIDArr = h5read(inpath, "cellIDArr")
    if ml==true
        cellIDArr = cellIDArr[25:275, 25:275, 5:160]
        # cellStateArr = cellStateArr[25:275, 25:275, 5:160]
        dx = 0.002
    else
        cellStateArr = h5read(inpath, "cellStateArr")
        dx = h5read(inpath, "dx")
        cellIDArr = sLtrimCellID(cellIDArr, cellStateArr)
        cellStateArr = nothing
    end
    # grainArr = h5read(inpath, "grainArr")
    out = getConvInfo(cellIDArr, dx)
    varnames = Dict{}()
    varnames[@Name(ellipseLengths)] = out[1]
    varnames[@Name(grainRadius)] = out[2]
    varnames[@Name(ellipseAxis)] = out[3]
    varnames[@Name(countedVolume)] = out[4]
    varnames[@Name(ellipseVolume)] = out[5]
    varnames[@Name(AspectRatio)] = out[6]
    varnames[@Name(volumeWeightedAspectRatio)] = out[7]
    varnames[@Name(ellipseVolumeWeightedAspectRatio)] = out[8]
    # println(typeof(varnames))
    # println(length(varnames))
    # println(varnames)
    for key in keys(varnames)
        println(key)
        h5write(outpath, key, varnames[key])
    end
end
end

function runMultiLayerAnalysis(namefragment="")
    indir = mlindir
    outdir = mloutdir
    stats = Dict{String, Any}()
    for i in searchdir(indir,namefragment)
        inpath = joinpath(indir, i)
        simname = (replace(i, ".h5"=>""))*"_ellipses_sliced.h5"
        outpath = joinpath(outdir, simname)
        println("Saving to :", outpath)
        isfile(outpath) ? continue : nothing
        # h5write(outpath, "SimulationName", simname)
        cellIDArr = h5read(inpath, "cellIDArr")
        cellIDArr = cellIDArr[25:275, 25:275, 5:160]
        dx = 0.002
        stats[i] = getSlicedStats(cellIDArr, dx, i)

        # for key in keys(varnames)
        #     println(key)
        #     h5write(outpath, key, varnames[key])
        # end
    end
    return stats
end


function getConvInfo(cellIDArr, dx)
    eL, eDir, vol = getEllipsoid(cellIDArr, dx)
    print(size(eL), size(eDir), size(vol))
    cellIDArr = nothing
    eL[eL.>10^3 .|| eL .< -10^3] .= dx
    println(count(eL[:,1] .== -2))
    # vol[eL[:,1] .== -2] .= 0
    # eL[eL[:,1] .== -2, :] .= 0 # For all grains with only 1 cell
    eD = [norm(eL[i,:]) for i in 1:length(eL[:,1])]
    eVol = 4/3*pi.*[prod(eL[i,:]) for i in 1:length(eL[:,1])]
    eVol[eVol .== NaN .|| eVol .== Inf .|| eVol .== -1] .= 0
    shapeDev = eVol./vol
    AR = getAspectRatios(eL, dx)
    println(typeof(vol), typeof(eVol))
    # println(AR)
    # redVol = vol[vol .!= 0]
    # totVol = sum(redVol[redVol .< percentile(redVol, 99) .&& redVol .> percentile(redVol, 1)])
    totEVol = sum(eVol[eVol .< percentile(eVol, 99) .&& eVol .> percentile(eVol, 1)])
    println(percentile(AR[AR .!= 0],99))
    println(percentile(AR[AR .!= 0],1))
    println(percentile(vol[vol.!= 0],99))
    println(percentile(vol[vol.!= 0],1))
    arInds = findall(AR .!= 0 .&& AR .< percentile(AR[AR .!= 0],99) .&& AR .> percentile(AR[AR .!= 0],1) .&& vol .!= 0 .&& vol .< percentile(vol[vol.!= 0],99) .&& vol .> percentile(vol[vol.!= 0],1)) #Trimming outliers in 1st and 99th percentile, and removing null-entries
    println(length(arInds))
    vwAR = (AR[arInds].*vol[arInds])./sum(vol[arInds])
    earInds = findall(AR .!= 0 .&& AR .< percentile(AR[AR .!= 0],99) .&& AR .> percentile(AR[AR .!= 0],1) .&& eVol .!= 0 .&& eVol .< percentile(eVol[eVol.!= 0],99) .&& vol .> percentile(eVol[eVol.!= 0],1)) #Trimming outliers
    vwAREll = (AR[earInds].*eVol[earInds])./sum(eVol[earInds])
    return eL, eD, eDir, vol, eVol, AR, vwAR, vwAREll
end

function totTest()
    runConvAnalysis()
    checkConvCriteria()
end

function mlTest()
    runConvAnalysis("", true)
end

mlindir = "/net/nas/results/240523_sf_mlConvergence/inputs/"
mloutdir ="/net/nas/results/240523_sf_mlConvergence/outputs/"
slindir = "/net/nas/results/240523_sf_slConvergence/inputs/"
sloutdir = "/net/nas/results/240523_sf_slConvergence/outputs/"

mldict = Dict{}()
mldict["roi0_219w"] = "/net/nas/results/240516_sf_DebuggingCrystalGrowth/240516_1211_sf_Simulation/Results.h5"
mldict["roi0_285w"] = "/net/nas/results/240523_sf_DebuggingCrystalGrowth/240523_0637_sf_Simulation/Results.h5"
mldict["roi0_407w"] = "/net/nas/results/240518_sf_DebuggingCrystalGrowth/240518_1918_sf_Simulation/Results.h5"
mldict["roi2_219w"] = "/net/nas/results/240522_sf_DebuggingCrystalGrowth/240522_1729_sf_Simulation/Results.h5"
mldict["roi2_407w"] = "/net/nas/results/240523_sf_DebuggingCrystalGrowth/240523_0109_sf_Simulation/Results.h5"

# checkConvCriteria()
# checkConvCriteria(true)

# compareGrainSizes(true)

function mlGetData(stats, dataname, returnvalue=0, scenariogragment="", slicefragment="", exactslicename=false)
    if returnvalue == 1
        returnvector = Vector{Any}()
        namevector = Vector{String}()
    end
    for scenariokey in keys(stats)
        occursin(scenariogragment, scenariokey) ? nothing : continue
        scenario=stats[scenariokey]
        for slicekey in keys(scenario)
            occursin(slicefragment, slicekey) ? nothing : continue
            exactslicename && slicefragment .!== slicekey ? continue : nothing
            slice = scenario[slicekey]
            typeof(slice) == String ? continue : nothing
            datamatches = occursin.(dataname, keys(slice))
            count(datamatches .==1) == 0 ? println(keys(slice)) : nothing
            count(datamatches .==1) > 1 ? print(Tuple(keys(slice))[datamatches]) : nothing
            dataname = Tuple(keys(slice))[datamatches][1]
            # !(dataname in keys(slice)) ? println(keys(slice)) : nothing
            println("For scenario ", scenariokey, " and slice ", slicekey, "; ", string(dataname), ": ")
            returnvalue == 1 ? (append!(returnvector, [slice[dataname]]);push!(namevector, scenariokey*slicekey)) : nothing
            returnvalue == 0 ? print(slice[dataname]) : nothing
            #println("For scenario ", scenariokey, " and slice ", slicekey, "; meanGrainRadius:  ", mean(slice["grainRadius"]), " and medianGrainRadius: ", median(slice["grainRadius"]), " and meansqrtArea:", median(sqrt.(slice["countedArea"])))
        end
    end
    returnvalue == 1 ? (return returnvector, namevector) : nothing
end

function mlHits(stats, dataname, scenariofragment="", slicefragment="")
    mlGetData(stats, dataname, 2, scenariofragment, slicefragment)
end

# data407w = mlGetData("grainRadius", 1, "407", "ex050", true)

function plotGrainSizeVsHeight()
    # data = mlGetData("grainRadius", 1, "407", "ex050", true)
    data407 = mlGetData(stats, "grainRadius", 1, "407", "ez050", true)[1]
    data219 = mlGetData(stats, "grainRadius", 1, "219", "ez050", true)[1]
    plot219 = [mean(data219[1]), mean(data219[2])]
    plot407 = [mean(data407[1]), mean(data407[2])]
    xdir = [0, 2]
    a = plot(xdir, plot219)
    plot!(xdir, plot407)
    display(a)
end

function plotGrainSizeVsPower()
    dataroi0 = mlGetData(stats, "grainRadius", 1, "roi0", "ex025", true)
    println(dataroi0[2])
    dataroi0 = dataroi0[1]
    dataroi2 = mlGetData(stats, "grainRadius", 1, "roi2", "ex025", true)
    println(dataroi2[2])
    dataroi2 = dataroi2[1]
    plotroi0 = [mean(dataroi0[1]), mean(dataroi0[2]), mean(dataroi0[3])]
    plotroi0med = [median(dataroi0[1]), median(dataroi0[2]), median(dataroi0[3])]
    plotroi0stdp = [mean(dataroi0[1])+std(dataroi0[1]), mean(dataroi0[2])+std(dataroi0[2]), mean(dataroi0[3])+std(dataroi0[3])]
    plotroi0stdm = [mean(dataroi0[1])-std(dataroi0[1]), mean(dataroi0[2])-std(dataroi0[2]), mean(dataroi0[3])-std(dataroi0[3])]
    Proi0 = [219, 285, 407]
    a = plot(Proi0, plotroi0)
    plot!(Proi0, plotroi0med)
    # plot!(Proi0, plotroi0stdp)
    # plot!(Proi0, plotroi0stdm)
    plotroi2 = [mean(dataroi2[1]), mean(dataroi2[2])]
    plotroi2med = [median(dataroi2[1]), median(dataroi2[2])]
    # plotroi2stdp = [mean(dataroi2[1])+std(dataroi2[1]), mean(dataroi2[2])+std(dataroi2[2])]
    # plotroi2stdm = [mean(dataroi2[1])-std(dataroi2[1]), mean(dataroi2[2])-std(dataroi2[2])]
    Proi2 = [219, 407]
    plot!(Proi2, plotroi2)
    plot!(Proi2, plotroi2med)
    ylims!(0,1.2*maximum(maximum.([plotroi2,plotroi0])))
    # plot!(Proi2, plotroi2stdp)
    # plot!(Proi2, plotroi2stdm)
    display(a)
end

function plotSlicedGrainSize(low,high, inbd=true, snippet="")
    inbd ? xlen = range(low, high) .* 2 : xlen = range(50,250).*2
    plot_font = "Computer Modern"
    default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)
    a = plot(legend=:topright, xlabel=L"Build height in $\mu m$", ylabel=L"Mean grain area in $\mu m^2$")
    currmax = 0
    returnvecea = Vector{Any}()
    returnvecca = Vector{Any}()
    for name in searchdir(mlindir, snippet)
        cellIDArr = h5read(joinpath(mlindir, name), "cellIDArr")
        outputea, outputca = zeros(length(xlen)), zeros(length(xlen))

        # Threads.@threads for i = low:high
        #     el, eDir, area = getEllipse(cellIDArr[50:250,50:250,i], 0.002)
        #     eDir = nothing
        #     el = el[el .!= 0 .&& el .< 10^4]
        #     ea = pi.*[prod(el[i,:]) for i in 1:length(el[:,1])]
        #     area = area[area .!= 0 .&& area .< percentile(area, 98)]
        #     outputea[i-low+1] = mean(ea)
        #     outputca[i-low+1] = mean(area)
        # end

        if inbd
            Threads.@threads for i = low:high
                area = getEllipseSimple(cellIDArr[50:250,50:250,i], 0.002)[3]
                eDir = nothing
                # el = el[el .!= 0 .&& el .< 10^4]
                # ea = pi.*[prod(el[i,:]) for i in 1:length(el[:,1])]
                area = area[area .!= 0 .&& area .> percentile(area[area .!= 0], 50)]
                # outputea[i-low+1] = mean(ea)
                outputca[i-low+1] = mean(area)
            end
        else
            Threads.@threads for i = 50:250
            area = getEllipseSimple(cellIDArr[50:250,i,low:high], 0.002)[3]
            eDir = nothing
            # el = el[el .!= 0 .&& el .< 10^4]
            # ea = pi.*[prod(el[i,:]) for i in 1:length(el[:,1])]
            area = area[area .!= 0 .&& area .> percentile(area[area .!= 0], 50)]
            # outputea[i-low+1] = mean(ea)
            outputca[i-50+1] = mean(area)
            end
        end


        

        # Threads.@threads for i = 50:250
        #     el, eDir, area = getEllipse(cellIDArr[50:250,i,low:high], 0.002)
        #     eDir = nothing
        #     el = el[el .!= 0 .&& el .< 10^4]
        #     ea = pi.*[prod(el[i,:]) for i in 1:length(el[:,1])]
        #     area = area[area .!= 0 .&& area .< percentile(area, 98)]
        #     outputea[i-50+1] = mean(ea)
        #     outputca[i-50+1] = mean(area)
        # end

        # output = [mean(sqrt.(prod.(eachrow(clamp.(getEllipse(cellIDArr[50:250,50:250,i], 0.002)[1],0.001,10000))))) for i=low:high]
        # output = [mean(getEllipse(cellIDArr[50:250,i,low:high], 0.002)[3]) for i=50:250]s
        # output = [maximum(getEllipse(cellIDArr[50:250,50:250,i], 0.002)[3]) for i=low:high]
        # out = [getEllipse(cellIDArr[50:250,i,low:high], 0.002)[3] for i=50:250]
        # out = [getEllipse(cellIDArr[50:250,50:250,i], 0.002)[3] for i=low:high]
        # output = [mean(layer[layer .> 100]) for layer in out] #threshhold of grain size 50 to match Weiyu results
        # plot!(xlen, outputea, label=name*" ellipse")
        plot!(xlen, outputca, label=name*" counted")
        print(name)
        currmax = max(currmax, maximum(outputea))
        append!(returnvecea, [outputea])
        append!(returnvecca, [outputca])
    end
    ylims!(0,1.2*currmax)
    display(a)
    return returnvecea, returnvecca
end

function plotConv()
    conv = checkConvCriteria()
    t219 = [conv[1], conv[7],conv[4]]
    t285 = [conv[2], conv[8],conv[5]]
    t407 = [conv[3], conv[9],conv[6]]
    res = [10, 5, 2]
    plot_font = "Computer Modern"
    default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)
    convPlot = plot(res,t219, yscale=:log10, legend=:topright, xlabel=L"Mesh Distance in $μm$", ylabel=L"Volume-Weighted Aspect Ratio $ \ $", label="219 W")
    plot!(res,t285, label="285 W")
    plot!(res,t407, label="407 W")
    xaxis!(xflip=true)
    # ylims!(0, 20)
    return convPlot
end

function plotPanzerFrac()
    plot_font = "Computer Modern"
    default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)
    # convPlot = plot(res,t219, yscale=:log10, legend=:topright, xlabel=L"Mesh Distance in $μm$", ylabel=L"Volume-Weighted Aspect Ratio $ \ $", label="219 W")
    # xdata = [267, 290, 296, 308, 333, 377, 400]
    # ydata = [2.62, 2.6, 2.95, 2.52, 2.4, 2.2, 2.1] #Trim 1 outlier
    xdata = [412, 390, 340, 310, 265, 295, 295] 
    xerrors = [85, 50, 60, 10, 5, 45, 25]
    ydata = [2.08, 2.2, 2.4, 2.5, 2.62, 2.6, 2.95]
    yerrors = [0.3, 0.25, 0.22, 0.25, 0.45, 0.4, 0.35]
    xlow = xdata .- xerrors
    xhigh = xdata .+ xerrors
    ylow = ydata .- yerrors
    yhigh = ydata .+ yerrors
    totalx = vcat(vcat(xdata,xlow), xhigh)
    totaly = vcat(vcat(ydata,ylow), yhigh)
    scatterplot = scatter(xdata, ydata, xlabel=L"Mean grain size $A_g$ in $μm^2$", ylabel=L"Crack length $L_c$ in $mm$", label="Experimental Data")
    # xdata = totalx
    # ydata = totaly
    # xdata = [267, 290, 308, 333, 377, 400] #1st try, based on density
    # ydata = [2.62, 2.6, 2.52, 2.4, 2.2, 2.1] #1st try
    # plot!(xdata, ydata)
    # plot!(xdata, ydata, xerr=[10,10,20,20,20, 10], yerr=[0.5, 0.3, 0.2, 0.5, 0.1, 0.2]) # Need to get error from the plots, these were just test values
    @. model(x,p) = p[1] + p[2]*x
    fit = curve_fit(model, xdata, ydata, [3, -0.001])
    xrange = range(minimum(totalx),maximum(totalx), 100)
    trendline = fit.param[1] .+fit.param[2].*xrange
    # combinedplot = scatterplot
    # plot!(xrange, trendline, label = "Best fit (linear): "*(@sprintf "%.3f" fit.param[2])*"x + "*(@sprintf "%.3f" fit.param[1]))
    # ylims!(0, 1.1*maximum(totaly))
    display(scatterplot)
    return scatterplot, fit.param
    # return [scatterplot, combinedplot]
    # line = model(xdata, )
end

function plotWeiyu()
    plot_font = "Computer Modern"
    default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)
    roi0_weiyu = [340, 430, 480]
    roi1_weiyu = [390, 495, 570]
    roi2_weiyu = [500, 670, 720]

    # roi0=[218.15, 203.2] #plotSlicedGrainSize(20,10,50,"3.5"), in yz plane
    # roi2=[176.48, 185]
    
    roi0 = [138.46, 163.88]
    roi2 = [157.62, 170.85] # plotSlicedGrainSize(20,10,50,"3.5"), in xy plane
    mypowers = [219, 407]
    powers = [219, 285, 407]
    explot = plot(powers,roi0_weiyu, legend=:bottomright, color=:orange, linestyle=:dash, xlabel=L"Laser Power in W $ \ $", ylabel=L"Mean grain area in $ \mu m^2 $", label="Scenario 0 - Experimental")
    plot!(powers, roi1_weiyu, color=:green,  linestyle=:dash, label="Scenario 1 - Experimental")
    plot!(powers, roi2_weiyu, color=:blue,  linestyle=:dash, label="Scenario 2 - Experimental")
    scatter!(powers, roi0_weiyu, color=:black)
    scatter!(powers, roi1_weiyu, color=:black)
    scatter!(powers, roi2_weiyu, color=:black)

    simplot = twinx()
    plot!(simplot, mypowers, roi0,color=:orange,legend=:topleft ,label="Scenario 1 - Simulative")
    plot!(simplot, mypowers, roi2,color=:blue,label="Scenario 2 - Simulative")
    scatter!(simplot, mypowers, roi0, color=:black)
    scatter!(simplot, mypowers, roi2, color=:black)

    ylims!(explot, 0.7*minimum(roi0_weiyu), 1.2*maximum(roi2_weiyu))
    ylims!(simplot, 0.7*minimum(roi0), 1.2*maximum(roi2))
    display(explot)
    # newy = twinx()
    # plot!(newy, mypowers, roi0)

    # roiplot = plot(mypowers, roi0,color=:red,legend=:bottomright, xlabel=L"Laser Power in W $ \ $", ylabel=L"Mean grain area in $ \mu m^2 $",label="Scenario 1 - Simulative")
    # plot!(mypowers, roi2,color=:green,label="Scenario 2 - Simulative")
    # scatter!(mypowers, roi0, color=:black)
    # scatter!(mypowers, roi2, color=:black)
    # ylims!(0.7*minimum(roi0), 1.2*maximum(roi2))
    # roi0 = []
    # roi2 = 
    return explot
end

function getVolume(snippet="")
    means = Vector{}()
    for filename in searchdir(mlindir,snippet)
        cellIDArr = h5read(joinpath(mlindir, filename), "cellIDArr")[50:250,50:250,25:100]
        grains = unique(cellIDArr[cellIDArr .!= 0])
        gsize = zeros(Int, length(grains))
        Threads.@threads for g in eachindex(grains)
            gsize[g] = count(cellIDArr .== grains[g])
            # push!(gsize, 
        end
        push!(means, [gsize[gsize .!= 0].*(0.002*1000)^3])
    end
    return means
    end

