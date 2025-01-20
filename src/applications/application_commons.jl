#some common structs and functions for all applications
using Lux

# Create a struct to hold all the constant parameters for the SWE2D RHS function
# These parameters are "extra", meaning, they are not part of the inversion/sensitivity analysis parameters, 
# such as Manning's n, zb, inletQ, etc. (These are part of the params_array)
struct SWE2D_Extra_Parameters{T<:Real}
    case_path::String
    active_param_name::String
    settings::Hydrograd.ControlSettings
    bInPlaceODE::Bool
    my_mesh_2D::Hydrograd.mesh_2D
    nodeCoordinates::Matrix{T}
    srh_all_Dict::Dict
    boundary_conditions::Hydrograd.BoundaryConditions2D
    swe_2D_constants::Hydrograd.swe_2D_consts
    ManningN_cells::Vector{T}
    inletQ_Length::Vector{Vector{T}}
    inletQ_TotalQ::Vector{T}
    exitH_WSE::Vector{T}
    wse::Vector{T}
    wstill::Vector{T}
    h::Vector{T}
    hstill::Vector{T}
    wse_ghostCells::Vector{T}
    wstill_ghostCells::Vector{T}
    h_ghostCells::Vector{T}
    hstill_ghostCells::Vector{T}
    zb_cells::Vector{T}
    zb_ghostCells::Vector{T}
    zb_faces::Vector{T}
    S0_cells::Matrix{T}
    S0_faces::Vector{Vector{Vector{T}}}
    b_dry_wet::BitVector
    b_Adjacent_to_dry_land::BitVector
    b_Adjacent_to_high_dry_land::BitVector
    ude_model::Union{Nothing, Lux.Chain}
    ude_model_params::Union{Nothing, NamedTuple}
    ude_model_state::Union{Nothing, NamedTuple}
end
