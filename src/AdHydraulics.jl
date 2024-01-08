module AdHydraulics

include("functions.jl")
include("mesh/mesh_SRHGEOM.jl")
include("mesh/mesh_GMSH2D.jl")

export greet_your_package_name
export readSRHGEOM
export readGmshMesh


end
