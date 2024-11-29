module AdHydraulics

include("functions.jl")

#swe_1D
include("mesh/mesh_1D.jl")
include("constants/swe_1D_constants.jl")
include("solver/swe_1D_solvers.jl")
include("util/swe_1D_tools.jl")

#swe_2D
include("mesh/mesh_SRHGEOM.jl")
include("mesh/mesh_GMSH2D.jl")

include("mesh/SRH_2D_Data.jl")

export greet_your_package_name

#swe_1D
export swe_1D_const
export mesh_1D
export initialize_mesh_1D

export Riemann_1D_hll!

export swe_1D_calc_total_water_volume
export swe_1D_save_results
export swe_1D_make_plots
export swe_1D_make_animation

#swe_2D
export read_srhgeom
export check_cell_nodes_counter_clockwise_srhgeom
export compute_mesh_properties_srhgeom
export read_srhmat

export read_gmsh_v2
export compute_mesh_properties

export SRH2D


end
