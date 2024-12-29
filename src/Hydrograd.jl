module Hydrograd

using PyCall
using JSON3
using JLD2
using LinearAlgebra
using Statistics 
using CSV
using DataFrames
using Printf
using Plots

#SciML related 
using SciMLSensitivity
using OrdinaryDiffEq

#Optimizers
using Optimization
using OptimizationOptimisers
#using OptimizationPolyalgorithms
#using OptimizationNLopt
#using Optim
#using OptimizationFlux
#using Flux

#AD engines
using Zygote
using ForwardDiff
using ReverseDiff
using Enzyme

#for Bayesian estimation
#using Turing


function print_banner()
    banner = """
    ╔════════════════════════════════════════════════════════════════════════╗
    ║                               Hydrograd                                ║
    ║              A Computational Hydrodynamics Package with                ║
    ║                        Automatic Differentiation                       ║
    ║------------------------------------------------------------------------║
    ║ Version    : 0.1.0                                                     ║
    ║ Author     : Xiaofeng Liu, PhD, P.E.                                   ║
    ║ Institution: Penn State University                                     ║
    ║ Email      : xzl123@psu.edu                                            ║
    ║ License    : MIT                                                       ║
    ║ Repository : https://github.com/psu-efd/Hydrograd.jl                   ║
    ╚════════════════════════════════════════════════════════════════════════╝
    """
    println(banner)
end

function greet_your_package_name()  
    return "Hello from Hydrograd!"
end

#applications
include("applications/inversion/process_inversion_results_2D.jl")

#constants
include("constants/swe_1D_constants.jl")
include("constants/swe_2D_constants.jl")

#controls
include("controls/control_settings_2D.jl")

#fvm
include("fvm/discretization/semi_discretize_2D.jl")
include("fvm/discretization/fvm_schemes_2D.jl")
include("fvm/discretization/Riemman_solvers/swe_1D_solvers.jl")
include("fvm/discretization/Riemman_solvers/swe_2D_solvers.jl")
include("fvm/boundary_conditions/bc_2D.jl")
include("fvm/initial_conditions/process_ICs_2D.jl")

#meshes
include("meshes/mesh_1D.jl")
include("meshes/mesh_2D.jl")

#ode solvers
include("ode_solvers/custom_ODE_solvers.jl")

#parameters
include("parameters/process_model_parameters_2D.jl")
include("parameters/process_ManningN_2D.jl")
include("parameters/process_bed_2D.jl")

#utilities
include("utilities/misc_tools.jl")
include("utilities/process_SRH_2D_input.jl")
include("utilities/smooth_functions.jl")
include("utilities/swe_1D_tools.jl")
include("utilities/swe_2D_tools.jl")


export print_banner
export greet_your_package_name

# Applications
export process_inversion_results_2D

# Constants
export swe_1D_const
export swe_2D_const
export update_swe_2D_constants!

# Controls
export TimeSettings
export ForwardSimulationSettings
export InversionSettings
export parse_control_file

# FVM
#     boundary conditions
export BoundaryConditions2D
export initialize_boundary_conditions_2D
export process_all_boundaries_2d

#     discretization - Riemann solvers
export Riemann_1D_hll!
export Riemann_2D_Roe

#     discretization - fvm schemes
export update_ghost_cells_scalar
export compute_scalar_gradients
export compute_cell_gradient
export cells_to_faces_scalar

#     discretization - semi-discretization
export semi_discretize_2D

#     discretization - initial conditions
export setup_initial_condition!
export setup_ghost_cells_initial_condition!

# Meshes
export mesh_1D
export mesh_2D
export initialize_mesh_1D
export initialize_mesh_2D

# ODE Solvers
export custom_ODE_solve
export custom_ODE_update_cells

# Parameters
export process_model_parameters_2D
export process_ManningN_2D
export process_bed_2D

# Utilities
export update_1d_array
export process_SRH_2D_input 
export smooth_abs
export smooth_sign
export smooth_max
export smooth_min
export smooth_sqrt
export smooth_pow
export swe_1D_calc_total_water_volume
export swe_1D_save_results
export swe_1D_make_plots
export swe_1D_make_animation
export swe_2D_calc_total_water_volume
export swe_2D_save_results_SciML
export swe_2D_save_results_custom

export export_to_vtk_2D


end
