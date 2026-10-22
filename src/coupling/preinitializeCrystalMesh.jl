using AMMicrostructure: updateInterfaceCells
using Random

function preinitializeCrystalMesh(p)
    nx, ny, nz, numGrains = p.nx, p.ny, p.nz, p.numGrains

    cellIDArr = zeros(Int32, (nx,ny,nz)) # Index of associated grain in grainArr
    grainArr = zeros(Float32, (6, nx*ny*nz)) #coords x,y,z : angles phi1,phi2,phi3
    print("would start voronoiInit")
    cellIDArr, grainArr = voronoiInit(grainArr, cellIDArr, p)

return cellIDArr, grainArr
end

function voronoiInit(grainArr, cellIDArr, p)
    if p.multilayer == true
        lengthx, lengthy= p.lengthx, p.lengthy
        lengthz = 5*p.dx
        nx, ny = p.nx, p.ny
        nz = 5
    else
        lengthx, lengthy, lengthz = p.lengthx, p.lengthy, p.lengthz
        nx, ny, nz = p.nx, p.ny, p.nz
    end
    print(nx, " ", ny, " ", nz, "\n")
    avgGrainVolume = (0.01)^3 #Expressed in mm just like all sizes.
    avgGrainSize = floor(Int, avgGrainVolume/p.dx^3)
    numGrains = floor(Int, lengthx*lengthy*lengthz/(avgGrainVolume))
    grainArr[1,1:numGrains] = rand(Float64, numGrains)*lengthx
    grainArr[2,1:numGrains] = rand(Float64, numGrains)*lengthy
    grainArr[3,1:numGrains] = rand(Float64, numGrains)*lengthz

    grainArr[4:6,1:numGrains] = rand(Float64, (3, numGrains)).*90 .-44
    hBIni = floor(Int,cbrt(avgGrainSize))
    println("HBini: ", hBIni)
    counters = zeros(Int, Threads.nthreads())
    Threads.@threads for grain in 1:numGrains
            hB = rand(1:hBIni)
            i,j,k = floor(Int, grainArr[1,grain]/p.dx)+1, floor(Int, grainArr[2,grain]/p.dx)+1, floor(Int, grainArr[3,grain]/p.dx)+1
            cellIDArr[(max((i-hB),1)):(min((i+hB),nx)), (max((j-hB),1)):(min((j+hB),ny)), (max((k-hB),1)):(min((k+hB),nz))] .= grain
            counters[Threads.threadid()] += 1
            if sum(counters) % floor(Int, numGrains/100) == 0
                print("Progress: ", floor(Int, sum(counters)/numGrains*100), "%\n")
            end
        end
    grainArr = grainArr[4:6,:]
    cellStateArr = zeros(Int8, nx, ny, nz)
    cellStateArr[cellIDArr[:,:,1:nz] .== 0] .= 2
    b = count(cellStateArr .== 2)
    while b > 0
        updates = updateInterfaceCells(cellStateArr, view(cellIDArr, :,:,1:nz), 0.7, nx, ny, nz)
        b = count(cellStateArr[:,:,:] .== 2)
        print("Progress: ", floor(Int, (1-b/(nx*ny*nz))*100), "%\n")
    end
    return cellIDArr, grainArr

end
