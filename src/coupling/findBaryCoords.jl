using LinearAlgebra
using HDF5
using AMMicrostructure: load_array, wait_for_key, to_dict
using Printf

function checkIfInsideElem(nodes, testPoint) 
    A = nodes[:, 1]
    B = nodes[:, 2]
    C = nodes[:, 3]
    D = nodes[:, 4]
    Vtot = det([B-A C-A D-A])
    V1 = det([testPoint-A C-A D-A])
    V2 = det([B-A testPoint-A D-A])
    V3 = det([B-A C-A testPoint-A])
    baryCoords = [0, V1, V2, V3]/Vtot
    baryCoords[1] = 1 - sum(baryCoords[2:4])
    any(baryCoords .< 0) || any(baryCoords .>1) ? (return (0, [0, 0, 0, 0])) : (return (1, baryCoords))
end

function calcBaryCoords(baryCoordsArrayPtr, elemArrayPtr, nodesPtr, numNodes, elemsPtr, konLength, ipkonPtr, ipkonLength, p) #C facing function with baryCoordsarray and elemArray changed
    print(typeof(baryCoordsArrayPtr), " ", typeof(elemArrayPtr), " ", typeof(nodesPtr), " ", typeof(elemsPtr), " ", typeof(ipkonPtr), "\n")

    resDir = "/net/nas/results/240514_sf_CrystalMeshInterpolation/multiLayer/roi0/"
    resName = joinpath(resDir,p.filename*"_geom.h5")
    h5write(resName, "test", 0)
    nx, ny, nz, dx = p.nx, p.ny, p.nz, p.dx
    baryCoordsArray = zeros(Float32, nx*ny*nz*4)
    load_array(baryCoordsArrayPtr, nx*ny*nz*4, baryCoordsArray)

    elemArray = zeros(Int32, nx*ny*nz*4)
    load_array(elemArrayPtr,nx*ny*nz,elemArray)

    nodes = zeros(Float64, 3*numNodes)
    load_array(nodesPtr, 3*numNodes, nodes)

    elems = zeros(Int, konLength)
    load_array(elemsPtr, konLength, elems)

    ipkon = zeros(Int, ipkonLength) #Starting location of each element in the elems array
    load_array(ipkonPtr, ipkonLength, ipkon)

    baryCoordsArray = reshape(baryCoordsArray, (nx,ny,nz,4))
    elemArray = reshape(elemArray,(nx,ny,nz, 4))
    nodes = reshape(nodes, (3, :))
    inWhichElem = zeros(Int, nx, ny, nz)

    Xmin, Ymin, Zmin = minimum(nodes, dims=2)

    numElems = ipkonLength
    print("NumElems: ", numElems, "NumNodes: ", numNodes, "nx: ", nx, "ny: ", ny, "nz: ", nz, "dx: ", dx, "Xmin: ", Xmin, "Ymin: ", Ymin, "Zmin: ", Zmin)
    wait_for_key("Press any key to continue")
    count = 0    

    Threads.@threads for elem in 1:numElems #export JULIA_NUM_THREADS=8 is configured in ~/.bashrc and julia_interface.c -> TODO: Check if it needs to be in bashrc
        ipkon[elem] < 0 ? continue : nothing
        currElemNodes = nodes[:, (elems[ipkon[elem]+1:ipkon[elem]+4])] #Note indexing here, +1 to +4 as first element is metadata
        elem_nx_min, elem_ny_min, elem_nz_min = max.(floor.(Int, (minimum(currElemNodes, dims=2) .- [Xmin, Ymin, Zmin])./dx .-1), 1)
        elem_nx_max, elem_ny_max, elem_nz_max = min.(ceil.(Int, (maximum(currElemNodes, dims=2) .- [Xmin, Ymin, Zmin])./dx .+1), [nx, ny, nz]) 
        for k in elem_nz_min:elem_nz_max
            for j in elem_ny_min:elem_ny_max
                for i in elem_nx_min:elem_nx_max
                    elemArray[i,j,k,1] != 0 ? continue : nothing # Test 
                    currPoint = [(i)/(nx+2)*p.lengthx+Xmin, (j)/(ny+2)*p.lengthy+Ymin, (k)/(nz+2)*p.lengthz+Zmin]
                    insideElem, baryCoords = checkIfInsideElem(currElemNodes, currPoint)
                    if(insideElem == 1)
                        round(sum(baryCoords),digits=1) != 1.0 ? print(i,j,k) : nothing # Test if baryCoords are valid (sum to 1)
                        elemArray[i,j,k, :] = elems[ipkon[elem]+1:ipkon[elem]+4]
                        baryCoordsArray[i,j,k,:] = baryCoords
                        inWhichElem[i,j,k] = elem
                    end
        end end end
        count += 1
        mod(count,100)==0 ? print("Percentage complete:" , count/numElems*100, "\n") : nothing
    end

    print("Number of elements without a point allocated: ", sum(elemArray .== 0), "\n")

    h5write(resName, "baryCoordsArray", baryCoordsArray)
    h5write(resName, "elemArray", elemArray)
    h5write(resName, "config", @sprintf("%s",(to_dict(p))))
    h5write(resName, "originCoords", [Xmin, Ymin, Zmin])

    return convert.(Float64, baryCoordsArray), convert.(Int, elemArray)
end

