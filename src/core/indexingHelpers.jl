# Indices for 1st neumann neighborhood
N1Ind = [
    0 0 1; 0 0 -1;
    0 1 0; 0 -1 0;
    1 0 0; -1 0 0
]

#To generate a neighborhood:
function getNeighborhood(size, type)
    inds = [CartesianIndex(i,j,k) for i=-size:size, j=-size:size, k=-size:size]
    inds = vec(reshape(inds, 1, :))
    c = Vector{CartesianIndex{3}}()
    for i in 1:length(inds)
        closestLeft = argmin([abs(inds[j][1])+abs(inds[j][2])+abs(inds[j][3]) for j in eachindex(inds)])
        append!(c,[inds[closestLeft]])
        deleteat!(inds,closestLeft)
    end
    type == 0 ? (return c[2:length(c)]) : (return c[2:((2*(size-1)+1)^3+6)])
end

function getMooreN(size)
    return getNeighborhood(size, 0)
end

function toCartesianInd(A::AbstractArray)
    indVec = Vector{CartesianIndex{3}}()
    for i=1:size(A,1)
        append!(indVec, [CartesianIndex(A[i,:]...)])
    end
    return indVec
end

# Octahedron face normals
OctoNormVec = (1/sqrt(3)*[-1 -1 -1; -1 -1 1; -1 1 -1; -1 1 1; 1 -1 -1; 1 -1 1; 1 1 -1; 1 1 1])

# Gets the octahedron coordinates of a cell, either retrieving or initializing them.
function cellCoords(i,j,k,p, cellOctaArr)
    all(cellOctaArr[1:3,i,j,k] .== 0) ? coords = [(i-0.5)*p.dx, (j-0.5)*p.dx, (k-0.5)*p.dx] : coords = cellOctaArr[1:3,i,j,k]
    return coords
end

# Rotates a vector by phi1, phi2, phi3 in degrees, in Z-X-Z order of rotation. To reverse rotation, negate the angles AND reverse their order
function rotVector(vector, phi1, phi2, phi3) 
    R = [-sind(phi1)*sind(phi3)*cosd(phi2)+cosd(phi1)*cosd(phi3) -sind(phi1)*cosd(phi3)-sind(phi3)*cosd(phi1)*cosd(phi2)  sind(phi2)*sind(phi3);
                   sind(phi1)*cosd(phi2)*cosd(phi3)+sind(phi3)*cosd(phi1) -sind(phi1)*sind(phi3)+cosd(phi1)*cosd(phi2)*cosd(phi3) -sind(phi2)*cosd(phi3);
                                                    sind(phi1)*sind(phi2)                                   sind(phi2)*cosd(phi1)             cosd(phi2)] #Euler rotation in Z-X-Z order of rotation
return R*vector
end
