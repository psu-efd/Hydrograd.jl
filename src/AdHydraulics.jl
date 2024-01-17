module AdHydraulics

include("functions.jl")

#swe_1D
include("mesh/mesh_1D.jl")
include("field/swe_1D_fields.jl")
include("solver/swe_1D_solvers.jl")

#swe_2D
include("mesh/mesh_SRHGEOM.jl")
include("mesh/mesh_GMSH2D.jl")

export greet_your_package_name

#swe_1D
export swe_1D_parameters
export Boundary_Type_Name, wall, zeroGradient, inletQ, exitH
export Boundary_1D
export mesh_1D
export swe_1D_fields
export initialize_mesh_1D
export initialize_swe_1D_fields

export Riemann_1D_hll
export swe_1D_rhs!

#swe_2D
export readSRHGEOM
export readGmshMesh


end
