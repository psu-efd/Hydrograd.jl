#1D mesh struct: only x coordinates, no bed elevation zb. Bed elevation information is external 
#to the mesh_1D struct so it can be a differentiable variable (not sure whether my understanding 
#is correct. It seems currently Zygote and Enzyme don't support struct with both differentiable
#and non-differentiable variables. This could change in the futre.)
mutable struct mesh_1D
    nCells::Int64     # Number of cells
    nFaces::Int64     # Number of faces (= nCells + 1)
    L::Float64  # Length of the domain
    xCells::Vector{Float64}  #x coordinates of cell centers 
    xFaces::Vector{Float64}  #x coordinates of faces (points in 1D) 

    dx::Float64  # Grid spacing (cell size)
end


# Function to initialize the mesh_1D struct
function initialize_mesh_1D(nCells, L) 
    dx = L / nCells
    nFaces = nCells + 1

    xCells = [ dx/2.0 + dx*i for i=0:(nCells-1) ]

    xFaces = [ dx*i for i=0:(nFaces-1) ]

    return mesh_1D(nCells, nFaces, L, xCells, xFaces, dx)
end

