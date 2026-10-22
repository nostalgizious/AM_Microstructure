# Cell-state updates, nucleation and decentered-octahedron growth.
using LinearAlgebra
using Random

include("indexingHelpers.jl")
include("configFunctions.jl")


# Pure julia function handling crystal growth
function updateCrystalMeshReduced!(cellIDArr, cellStateArr, cellOctaArr, grainArr, temps, prevTemps, logVars, params, nx, ny, currNz, currZmin, currZmax, currHeight, numGrains, dt, istep;rng=Random.default_rng())
    println(params)
    maxPermittedTempChange = (params.Tliquidus-params.Tsolidus-1)/30
    maxCooled = getMaxCooled(cellStateArr, temps, prevTemps, maxPermittedTempChange, params)
    subStepsTemp = ceil(Int, maxCooled/maxPermittedTempChange) # Max temp reduction of liquid / max allowed cooling in 1 step => substeps needed *2 for margin
    subTimeSteps = clamp(subStepsTemp, 4, 400) # 40 is max substeps, very high, but only needed for very high cooling rates which occur ~5-10% of the time
    println("Substeps set to:", subTimeSteps," Substeps needed for dendrite growth: ", 0, " Substeps needed for temperature change: ", subStepsTemp)

    timeUpdateCellStates, timeNucleation, timeUpdateInterfaceCells, timeUpdateGrowingCells, timeToUpdateDendriteLengths = 0, 0, 0, 0, 0
    numCellsNucleated, numInterfaceGrowth, numOctahedralGrowth = 0,0,0
    M2Ind = getMooreN(2)
    
    factor = 1/subTimeSteps
    tempDiffPerStep = similar(temps)
    Threads.@threads for i in eachindex(temps)
        tempDiffPerStep[i] = factor*(temps[i]-prevTemps[i])
    end
    temps = nothing

    Nmax=0
    faceCaptures = zeros(Int, 8)
    numCoolingTooQuickly = 0
    oldNumGrains = numGrains
    locTimestep = dt/subTimeSteps
    for timestep = 1:subTimeSteps
        timestep == 1 ? (locPrevTemps = prevTemps; prevTemps = nothing) : locPrevTemps = locTemps
        locTemps = tempTimeInterpolation(locPrevTemps, tempDiffPerStep)


        preTime = time()
        updateCellStates(cellStateArr, cellIDArr, cellOctaArr, locPrevTemps, locTemps, params, logVars, currNz)
        timeUpdateCellStates += time()-preTime
        print("Cell states updated")
        
        preTime = time()
        numGrains,currNmax = nucleateCells(cellStateArr, grainArr, cellIDArr, locTemps, locPrevTemps, params, numGrains, logVars;rng)
        Nmax = max(Nmax, currNmax)
        print("Nucleation done")
        timeNucleation += time()-preTime
        
        preTime = time()
        numInterfaceGrowth += updateInterfaceCells(cellStateArr, cellIDArr, 0.99, nx, ny, currNz;rng)
        timeUpdateInterfaceCells += time()-preTime
        print("Interface cells updated")

        preTime = time()
        updateDendriteLengths(cellStateArr, cellOctaArr, locTemps, params, locTimestep)
        timeToUpdateDendriteLengths += time()-preTime
        

        preTime = time()
        newOctahedralCells, newInterfaceCells, whichFacesCaptured = updateCellGrowth(cellStateArr, cellIDArr, cellOctaArr, grainArr, params, M2Ind, toCartesianInd(N1Ind), OctoNormVec, currNz;rng)
        numOctahedralGrowth += newOctahedralCells
        numInterfaceGrowth += newInterfaceCells
        faceCaptures += whichFacesCaptured

        timeUpdateGrowingCells += time()-preTime
        print("Growing cells updated")
    end
    temps = locTemps
    locTemps, locPrevTemps = nothing, nothing
    numCellsNucleated = numGrains - oldNumGrains

    preTime = time()
    mmmFace, mmpFace, mpmFace, mppFace, pmmFace, pmpFace, ppmFace, pppFace = faceCaptures
    logVars["mmmFace"] = mmmFace
    logVars["mmpFace"] = mmpFace
    logVars["mpmFace"] = mpmFace
    logVars["mppFace"] = mppFace
    logVars["pmmFace"] = pmmFace
    logVars["pmpFace"] = pmpFace
    logVars["ppmFace"] = ppmFace
    logVars["pppFace"] = pppFace
    logVars["substeps"] = subTimeSteps
    logVars["numActiveSolids"] = count(cellStateArr .== 1)
    logVars["numLiquidCells"] = count(cellStateArr .== 3)
    logVars["numCoolingTooQuickly"] = numCoolingTooQuickly
    logVars["numGrains"] = numGrains
    logVars["numInterfaceGrowth"] = numInterfaceGrowth
    logVars["numOctahedralGrowth"] = numOctahedralGrowth
    logVars["numCellsNucleated"] = numCellsNucleated
    logVars["currZmin"] = currZmin
    logVars["currZmax"] = currZmax
    logVars["currHeight"] = currHeight
    logVars["maxNucleation"] = Nmax
    logVars["timeUpdateCellStates"] = timeUpdateCellStates
    logVars["timeNucleation"] = timeNucleation
    logVars["timeUpdateInterfaceCells"] = timeUpdateInterfaceCells
    logVars["timeUpdateGrowingCells"] = timeUpdateGrowingCells
    logVars["timeToLogVars"] = time()-preTime
    return numGrains
end     

function getMaxCooled(cellStateArr, temps, prevTemps, maxPermittedTempChange, p)
    numThreads = Threads.nthreads()
    maxCoolArr = zeros(numThreads)
    Threads.@threads :static for worker in 1:numThreads
        for i in (fld((worker-1)*length(temps),numThreads)+1):fld(worker*length(temps),numThreads)
            if cellStateArr[i]==3 && temps[i]<p.Tliquidus-maxPermittedTempChange
                maxCoolArr[worker] = max(maxCoolArr[worker],prevTemps[i]-temps[i])
            end
        end
    end
    return max(maximum(maxCoolArr), maxPermittedTempChange)
end
function tempTimeInterpolation(locPrevTemps, tempDiffPerStep)
    locTemps = similar(locPrevTemps)
    Threads.@threads for i in eachindex(locPrevTemps)
        locTemps[i] = locPrevTemps[i] + tempDiffPerStep[i]
    end
    return locTemps
end

# Updates cell states based on temperature
function updateCellStates(cellStateArr, cellIDArr, cellOctaArr, prevTemps, temps, p, logVars, currNz)    
    Threads.@threads for i in eachindex(temps)
        if temps[i] >= p.Tliquidus # Hotter than Tliquidus
            cellStateArr[i] = 3 # Liquid cell
            cellIDArr[i] = 0
            cellOctaArr[:,CartesianIndices(temps)[i]] .= 0 
        elseif cellStateArr[i] != 1 && (temps[i] > p.Tsolidus || cellIDArr[i] == 0) # All cells except freshly nucleated/grown ones, and between Tsolidus and Tliquidus are mushy  
            cellStateArr[i] = 2 # Mushy cell
            cellIDArr[i] = 0 # No grain assigned anymore. Do not remove entry in grainArr as it may belong to another cell as well
            cellOctaArr[:,CartesianIndices(temps)[i]] .= 0
        else
            nothing
        end
    end

    
end

# Nucleates new grains in mushy cells with a certain probability given by nucleationIntegral
function nucleateCells(cellStateArr, grainArr, cellIDArr, temps, prevTemps, p, numGrains, logVars;rng=Random.default_rng())
    nx, ny, nz, dx, Tliquidus, dTc, dTsigma, N0 = p.nx, p.ny, p.nz, p.dx, p.Tliquidus, p.dTc, p.dTsigma, p.N0

    candidates = findall(==(2),cellStateArr)
    isempty(candidates) && return numGrains,0.0
    draws = rand(rng,length(candidates))
    rates = zeros(length(candidates))
    Threads.@threads :static for index in eachindex(candidates)
        cell = candidates[index]
        rates[index] = min(nucleationIntegral(Tliquidus-temps[cell],Tliquidus-prevTemps[cell],p),0.005)
    end
    Nmax = maximum(rates)
    # Assign IDs and orientations in cell order, independent of task scheduling.
    for index in eachindex(candidates)
        draws[index] > rates[index] && continue
        cell = candidates[index]
        numGrains += 1
        cellIDArr[cell],cellStateArr[cell] = numGrains,1
        grainArr[:,numGrains] .= 90 .* rand(rng,3) .- 45
    end

    print("Nucleation rate max: ", Nmax, " Number of grains: ", numGrains, "\n")
    return numGrains, Nmax
end

# Determine cell specific nucleation rate adapted from Eq. 5.22 in (Beuerlein, 2023), originating from (Gür&Pan, 2009)
function nucleationIntegral(dT, dTprev, p) #TODO: Need to double check for more changes from 2D to 3D
    N0, dTc, dTsigma, dx = p.N0, p.dTc, p.dTsigma, p.dx
    # Composite Simpson quadrature for the original Gaussian law. Bound the
    # tails at 12 sigma (omitted mass < 4e-33), then use at least 32 panels.
    # N0's historical normalization is intentionally retained (no 1/dTsigma).
    dT > dTprev || return 0.0
    lo = max(dTprev, dTc - 12dTsigma)
    hi = min(dT, dTc + 12dTsigma)
    hi > lo || return 0.0
    panels = max(32, 2ceil(Int, 16(hi-lo)/dTsigma))
    h = (hi-lo)/panels
    gaussian(x) = exp(-0.5*((x-dTc)/dTsigma)^2)
    total = gaussian(lo) + gaussian(hi)
    for i in 1:panels-1
        total += (isodd(i) ? 4 : 2)*gaussian(lo+i*h)
    end
    return N0/sqrt(2pi) * h/3 * total * dx^3

end

# Updates all mushy cells that havent nucleated to be part of any grain in its 1st neumann neighborhood
function updateInterfaceCells(cellStateArr, cellIDArr, threshhold, nx, ny, currNz;rng=Random.default_rng())
    count(cellStateArr .== 0) == 0 ? (return 0) : nothing
    globInds = findall(cellStateArr .== 2)
    counters = Ref(0)
    if length(globInds) > 0
        for globInd in globInds
            i,j,k = globInd[1], globInd[2], globInd[3]
            for dir in 1:6  # Iterates 1st neumann neighborhood in a not very efficient way, but randomly to prevent bias. Hmm Randomness probably not desired, prefer 00-1 direction
                x, y, z = [i,j,k] .+ N1Ind[dir,:]
                x < 1 || x > nx || y < 1 || y > ny || z < 1 || z > currNz ? continue : nothing
                cellStateArr[x,y,z] > 1 || rand(rng)<threshhold || cellIDArr[x,y,z] == 0 ? continue : nothing # Skip neighbor if it isnt a Solid cell. rand()>probabilityOfSkipping -> To give dendritic growth more space. Does not necessarily skip, as it may interface several solid cells
                cellIDArr[i,j,k] = cellIDArr[x,y,z]
                cellStateArr[i,j,k] = 1
                counters[] += 1
                break
            end
        end
    end

    println("Cells updated: ", counters[])
    return counters[]
end

function updateDendriteLengths(cellStateArr, cellOctaArr, temps, p, locTimestep)
    any(==(2),cellStateArr) || return nothing
    # First Moore neighbourhood: all 26 face/edge/corner neighbours, clipped at
    # domain boundaries. A separable 3×3×3 dilation avoids scanning 26 neighbours
    # for every bulk solid cell. State 1 excludes the dilation's centre itself.
    nx,ny,nz = size(cellStateArr)
    along_x = Array{Bool}(undef,nx,ny,nz)
    along_xy = similar(along_x)
    Threads.@threads :static for k in 1:nz
        for j in 1:ny,i in 1:nx
            along_x[i,j,k] = cellStateArr[i,j,k] == 2 ||
                (i > 1 && cellStateArr[i-1,j,k] == 2) ||
                (i < nx && cellStateArr[i+1,j,k] == 2)
        end
    end
    Threads.@threads :static for k in 1:nz
        for j in 1:ny,i in 1:nx
            along_xy[i,j,k] = along_x[i,j,k] ||
                (j > 1 && along_x[i,j-1,k]) || (j < ny && along_x[i,j+1,k])
        end
    end
    Threads.@threads :static for k in 1:nz
        for j in 1:ny,i in 1:nx
            if cellStateArr[i,j,k] == 1 && (along_xy[i,j,k] ||
                (k > 1 && along_xy[i,j,k-1]) || (k < nz && along_xy[i,j,k+1]))
                cellOctaArr[4,i,j,k] += dendriteTipVelocity(p.Tliquidus-temps[i,j,k],p)*locTimestep
            end
        end
    end

    return 0
end

function updateCellGrowth(cellStateArr, cellIDArr, cellOctaArr, grainArr, p,
                          M2Ind::Vector{CartesianIndex{3}}, N1Ind::Vector{CartesianIndex{3}},
                          OctoNormVec, currNz;rng=Random.default_rng())
    # Capture against the state at the beginning of this sweep. New cells become
    # active only after all candidates have been examined (no chain reaction).
    faceCaptures = zeros(Int, 8)
    targets = findall(==(2),cellStateArr)
    isempty(targets) && return 0,0,faceCaptures
    captured_faces = zeros(UInt8,length(targets))
    sweep_seed = rand(rng,UInt64)
    # Source grains do not change during this sweep. Reuse their rotated faces.
    normals_by_grain = [rotVector(OctoNormVec',grainArr[:,id]...) for id in 1:maximum(cellIDArr)]
    Threads.@threads :static for index in eachindex(targets)
        target = targets[index]
        # Explicit integer mixing through Xoshiro seeding; no randomized hash().
        cell_rng = Random.Xoshiro(sweep_seed ⊻ UInt64(LinearIndices(cellStateArr)[target]))
        targetCentre = cellCoords(Tuple(target)..., p, cellOctaArr)
        # Retain the original preference for the six face-adjacent cells,
        # then consider the rest of the full second Moore neighborhood.
        for offset in vcat(shuffle(cell_rng,M2Ind[1:6]), shuffle(cell_rng,M2Ind[7:end]))
            source = target + offset
            checkbounds(Bool, cellStateArr, source) || continue
            cellStateArr[source] == 1 || continue
            cellIDArr[source] > 0 || continue
            Lnu = cellOctaArr[4, source]
            Lnu > 0 || continue
            sourceCentre = cellCoords(Tuple(source)..., p, cellOctaArr)
            normals = normals_by_grain[cellIDArr[source]]
            # An octahedron is the INTERSECTION of its eight half-spaces.
            # The face with the largest signed projection bounds this point.
            distances = normals' * (targetCentre - sourceCentre)
            face = argmax(distances)
            distances[face] <= Lnu || continue
            engulfCell(cellIDArr, cellStateArr, grainArr, cellOctaArr, p,
                       Tuple(source)..., Tuple(target)..., Lnu,
                       normals[:, face], Vector(OctoNormVec[face, :]))
            captured_faces[index] = face
            break
        end
    end
    # No target becomes a source until all worker tasks have joined.
    for index in eachindex(targets)
        face = captured_faces[index]
        face == 0 && continue
        cellStateArr[targets[index]] = 1
        faceCaptures[face] += 1
    end
    return sum(faceCaptures), 0, faceCaptures
end


#Model based on Kurz et al (polynomial approximation if computationally necessary)
function dendriteTipVelocity(dT, p)
    # Values for steel, copied from (Beuerlein, 2023), need to find own values, and read from config
    D_L, Gamma, m_L, k_0, C_0 = p.D_L, p.Gamma, p.m_L, p.k_0, p.C_0
    return 0.5*(D_L / (5.51 * (pi^2) * ((-m_L * (1 - k_0))^1.5) * Gamma)) * ((dT.^2.5) / (C_0^1.5)) #TODO: check different scaling here as calibration parameters, 0.5 for 200 mesh
end

function engulfCell(cellIDArr, cellStateArr, grainArr, cellOctaArr, p, i,j,k,x,y,z, Lnu, N, Nunrot)
    capturingGrainAngles = grainArr[:,cellIDArr[i,j,k]]
    capturingCellCoords = cellCoords(i,j,k,p,cellOctaArr)
    engulfedCellCoords = cellCoords(x,y,z,p,cellOctaArr)
    newCmu, newLmu = calcCentreOfEngulfedCell(capturingGrainAngles, capturingCellCoords, engulfedCellCoords, Lnu, N, Nunrot, p) # Save new centre of the engulfed cell, which is no longer exactly at it's corresponding grid point in cartesian CoSy
    cellOctaArr[1:3, x,y,z] = newCmu
    cellOctaArr[4, x,y,z] = newLmu
    cellIDArr[x,y,z] = cellIDArr[i,j,k]
end

#Implementation of 3D decentered octahedron algorithm(Gandin&Rappaz, 1997)
function calcCentreOfEngulfedCell(PhiNu, Cnu, Cmu, Lnu, N, Nunrot, p) 
    phi1,phi2,phi3 = PhiNu[1], PhiNu[2], PhiNu[3]

    
    # (ii) Calculate A = projection of cell centre Cmu onto face <111>
    # print("  N: ", N)
    # print("  Cmu: ", Cmu)
    # print("  Cnu: ", Cnu)
    # print("  Lnu: ", Lnu)
    
    # TODO: Check if this projection is correct
    d = Cmu .- Cnu
    A_unshifted = Lnu*N + (d .- dot(d,N)*N)  # points from origin rather than Cnu
    A = A_unshifted+Cnu #Projection of vector to captured cell centre, onto plane defined by normal vector N and with signed distance Lnu from Cnu
    # println(A)
    # (iii) Determine S1, corner of face <111>, which is the closest to A
    Aunrot = rotVector(A_unshifted, -phi3, -phi2, -phi1) 
    # println("An:", Aunrot)
    # Nunrot = rotVector(N,-phi1, -phi2, -phi3)
    Corners = round.(Nunrot, RoundFromZero).*Lnu*sqrt(3) #Rounds to nearest corner in an unrotated CoSy -> transforms e.g. 
    # println("Corners: ", Corners)
    # println("Cornerdist: ", Corners .- Aunrot)
    S1Ind = argmin([sum(abs2, Aunrot - [axis == j ? Corners[j] : 0.0 for axis in 1:3]) for j in 1:3]) #Find closest corner to A
    # println("Closestcorner: ", S1Ind)
    S1, S2, S3 = zeros(3), zeros(3), zeros(3)
    S1[S1Ind] = Corners[S1Ind]

    S2[(S1Ind)%3+1] = Corners[(S1Ind)%3+1] #Abitrary which is closer, does not matter, only S1 has to be closest
    S3[(S1Ind+1)%3+1] = Corners[(S1Ind+1)%3+1]
    # println(S1, S2, S3)
    S1, S2, S3 = rotVector(S1, phi1, phi2, phi3) .+ Cnu, rotVector(S2, phi1, phi2, phi3) .+ Cnu, rotVector(S3, phi1, phi2, phi3) .+ Cnu #Transform back to original CoSy
    # println("S1: ",S1)
    # println("S2: ", S2)
    # println("S3: ", S3)
    # (iv) Lengths ISi and JSi from I and J (projections of A onto segment edges S1S2 and S1S3to the vertices of the octahedron
    # println("I:", dot(A.-S1, normalize!(S2 .- S1)))
    # println("J:", dot(A.-S1, normalize!(S3 .- S1)))
    # I = dot(A.-S1, normalize!(S2 .- S1))*normalize!(S2.-S1)
    # J = dot(A.-S1, normalize!(S3 .- S1))*normalize!(S3.-S1)
    IS1 = dot(A.-S1, normalize!(S2 .- S1))
    JS1 = dot(A.-S1, normalize!(S3 .- S1))
    IS2 = norm(S2 .- S1) - IS1
    JS3 = norm(S3 .- S1) - JS1
    # IS1, IS2 =  norm(I .-S1), norm(I .-S2)
    # JS1, JS3 =  norm(J .-S1), norm(J .-S3)
    # println("IS1: ", IS1)
    # println("IS2: ", IS2)
    # println("JS1: ", JS1)
    # println("JS3: ", JS3)
    # (v) Calculate initial size of octahedron 
    L12 = 0.5*(min(IS1, sqrt(3)*p.dx) + min(IS2, sqrt(3)*p.dx))
    L13 = 0.5*(min(JS1, sqrt(3)*p.dx) + min(JS3, sqrt(3)*p.dx))
    # println("L12: ", L12)
    # println("L13: ", L13)
    Lmu = sqrt(2/3)*max(L12,L13)
    # println("Lmu: ", Lmu)
    # (vi) Calculate position of growth centre Cμ 
    newCmu = S1 .+sqrt(3)*Lmu*normalize!(Cnu .-S1)
    # println("Cmu:", Cmu)
    # print("  newCmu: ", newCmu, "\n")
    return newCmu, Lmu
end
