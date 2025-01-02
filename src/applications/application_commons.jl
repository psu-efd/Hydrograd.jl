#some common structs and functions for all applications

# Create a struct to hold all the constant parameters for the SWE2D RHS function
# These parameters are "extra", meaning, they are not part of the inversion/sensitivity analysis parameters, 
# such as Manning's n, zb, inletQ, etc. (These are part of the params_array)
struct SWE2D_Extra_Parameters{T}
    active_param_name::String
    settings::Hydrograd.ControlSettings
    my_mesh_2D::Hydrograd.mesh_2D
    nodeCoordinates::Matrix{T}
    srh_all_Dict::Dict
    boundary_conditions::Hydrograd.BoundaryConditions2D
    swe_2D_constants::Hydrograd.swe_2D_consts
    ManningN_cells::Vector{T}
    inletQ_Length::Vector{Vector{T}}
    inletQ_TotalQ::Vector{T}
    exitH_WSE::Vector{T}
    zb_cells::Vector{T}
    zb_ghostCells::Vector{T}
    zb_faces::Vector{T}
    S0::Matrix{T}
end
