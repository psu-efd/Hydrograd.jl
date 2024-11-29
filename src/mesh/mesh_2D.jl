#2D mesh struct: only non-mutable variables such as node coordinates
mutable struct mesh_2D
    nCells::Int64     # Number of cells
    nFaces::Int64     # Number of faces 
    nNodes::Int64     # Number of nodes
    
    xCells::Vector{Float64}  #x coordinates of cell centers 
    yCells::Vector{Float64}  #y coordinates of cell centers 
    xFaces::Vector{Float64}  #x coordinates of face centers 
    yFaces::Vector{Float64}  #y coordinates of face centers
    xNodes::Vector{Float64}  #x coordinates of nodes
    yNodes::Vector{Float64}  #y coordinates of nodes
    
    cellArea::Vector{Float64}  #cell area
    faceLength::Vector{Float64}  #face length

    faceNormalX::Vector{Float64}  #x component of face normal vector
    faceNormalY::Vector{Float64}  #y component of face normal vector

end


# Function to initialize the mesh_2D struct
function initialize_mesh_2D(nCells, L) 
    dx = L / nCells
    nFaces = nCells + 1

    xCells = [ dx/2.0 + dx*i for i=0:(nCells-1) ]

    xFaces = [ dx*i for i=0:(nFaces-1) ]

    return mesh_1D(nCells, nFaces, L, xCells, xFaces, dx)
end

