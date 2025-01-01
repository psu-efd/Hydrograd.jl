using Revise

using Dates

using JSON3

using JLD2

using IterTools

using PyCall

using Profile

#includet(joinpath(dirname(dirname(dirname(@__FILE__))), "src", "Hydrograd.jl"))
#includet(joinpath(dirname(dirname(dirname(@__FILE__))), "src", "Hydrograd.jl"))
using Hydrograd

#SciML
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

#using LinearAlgebra

using InteractiveUtils

#using Random
#Random.seed!(1234)

#hack
#get the current directory
current_dir = pwd()
cd("C:\\Users\\xzl123\\research\\Hydrograd.jl\\examples\\SWE_2D")

# Add some warnings/info for gradient computation
SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true

#Timing 
start_time = now()  # Current date and time

#include the files to define the variables
include("define_variables_swe_2D.jl")

#print the banner
Hydrograd.print_banner()

#directory to read from/save to (the same directory as the current directory)
#save_path = dirname(@__FILE__)
case_path = pwd()

# Read and parse control file
println("Reading control file...")
control_file = joinpath(case_path, "run_control.json")
settings = Hydrograd.parse_control_file(control_file)

#define a swe_2D_constants object with some values from the control file
swe_2D_constants = Hydrograd.swe_2D_consts(t=settings.time_settings.tspan[1], dt=settings.time_settings.dt, 
    tStart=settings.time_settings.tspan[1], tEnd=settings.time_settings.tspan[2], tspan=(settings.time_settings.tspan[1], settings.time_settings.tspan[2]))

#read data from SRH-2D hydro, geom, and material files; it aslo create the 2D mesh.
srh_all_Dict = Hydrograd.process_SRH_2D_input(settings, case_path)

#update swe_2D_constants based on the SRH-2D data
if settings.time_settings.bUse_srhhydro_time_settings
    Hydrograd.update_swe_2D_constants!(swe_2D_constants, srh_all_Dict)
end

#Get the 2D mesh 
my_mesh_2D = srh_all_Dict["my_mesh_2D"]

# define solution variables at cells
wse = zeros(my_mesh_2D.numOfCells)          #free surface elevation at cells 
h = zeros(my_mesh_2D.numOfCells)            #water depth at cells 
q_x = zeros(my_mesh_2D.numOfCells)          #q_x=hu at cells 
q_y = zeros(my_mesh_2D.numOfCells)          #q_y=hv at cells 

# define solution variables at ghost cells
wse_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #free surface elevation at ghost cells 
h_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)            #water depth at ghost cells 
q_x_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #q_x=hu at ghost cells 
q_y_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #q_y=hv at ghost cells 

#mesh related variables which might be mutable (should not be in struct mesh_2D)
#Get nodeCoordinates: Float64 2D array [numOfNodes, 3]
nodeCoordinates = srh_all_Dict["srhgeom_obj"].nodeCoordinates

#  setup bed elevation 
#If performing inversion on zb, zero out the bed elevation in nodeCoordinates.
if settings.bPerform_Inversion && settings.inversion_settings.active_param_names == ["zb"]
    nodeCoordinates[:, 3] .= 0.0
end

#setup bed elevation: compute zb at cell centers (zb_cells) from nodes, 
#then interpolate zb from cell to face (zb_faces) and ghost cells (zb_ghostCells),
#and compute the bed slope at cells (S0). 
#Note: If performing inversion on zb_cells, these values are not used in the inversion process. 
#      Instead, they are updated in the inversion process.
zb_cells, zb_ghostCells, zb_faces, S0 = Hydrograd.setup_bed(settings, my_mesh_2D, nodeCoordinates, case_path, true)

#define the true bed elevation at cells
zb_cells_truth = zeros(size(zb_cells))

#If performing forward simulation, make a copy of the bathymetry truth: zb_cell_truth; otherwise for 
#inversion and sensitivity analysis, it is zero (its value will be loaded from the forward simulation result in the inversion).
if settings.bPerform_Forward_Simulation
    zb_cells_truth = deepcopy(zb_cells)
end

#get the true Manning's n and inlet discharges
srhhydro_ManningsN_Dict = srh_all_Dict["srhhydro_ManningsN"]

#define the true Manning's n values
ManningN_values_truth = zeros(length(srhhydro_ManningsN_Dict))

#If performing forward simulation, make a copy of the Manning's n truth: ManningN_values_truth; otherwise for 
#inversion and sensitivity analysis, it is zero (its value will be loaded from the forward simulation result in the inversion).
if settings.bPerform_Forward_Simulation
    ManningN_values_truth = Float64[srhhydro_ManningsN_Dict[i] for i in 0:(length(srhhydro_ManningsN_Dict)-1)]
end

#get the true inlet discharges (could be nothing if no inletQ_BCs)
srhhydro_inletQ_Dict = srh_all_Dict["srhhydro_IQParams"]

#define the true inlet discharges
inlet_discharges_truth = nothing

if length(srhhydro_inletQ_Dict) > 0
    inlet_discharges_truth = [parse(Float64, srhhydro_inletQ_Dict[i][1]) for i in 1:(length(srhhydro_inletQ_Dict))]
end

#print the true values. The true parameter values are computed regardless of whether performing forward simulation or inversion
if settings.bVerbose
    println("True zb_cells = ", zb_cells_truth)
    println("True Manning's n values = ", ManningN_values_truth)
    println("True inlet discharges = ", inlet_discharges_truth)
end

#preprocess: create a 1D array to combine all model parameters for 2D shallow water equations
#params_vector: the 1D array of the active parameter (zb_cells_truth, ManningN_values_truth, or inlet_discharges_truth)
#active_param_name: the name of the active parameter (zb, ManningN, or Q)
params_vector, active_param_name = Hydrograd.setup_model_parameters_2D(settings, my_mesh_2D, srh_all_Dict, zb_cells_truth, ManningN_values_truth, inlet_discharges_truth)

#Initial setup of Manning's n at cells and ghost cells using the SRH-2D data (if performing inversion on Manning's n, ManningN_cells will be updated later in the inversion process)
ManningN_cells, ManningN_ghostCells = Hydrograd.setup_ManningN(settings, my_mesh_2D, srh_all_Dict)

#setup initial condition for wse, h, q_x, q_y at cells
Hydrograd.setup_initial_condition!(settings, my_mesh_2D, nodeCoordinates, wse, zb_cells, h, q_x, q_y, swe_2D_constants, case_path, true)

#setup initial condition for wse, h, q_x, q_y at ghost cells
Hydrograd.setup_ghost_cells_initial_condition!(settings, my_mesh_2D, wse, h, q_x, q_y, wse_ghostCells, h_ghostCells, q_x_ghostCells, q_y_ghostCells)

#create and preprocess boundary conditions: boundary_conditions only contains the static information of the boundaries.
boundary_conditions, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length, inletQ_TotalA, inletQ_DryWet,
exitH_WSE, exitH_H, exitH_A, wall_H, wall_A, symm_H, symm_A = Hydrograd.initialize_boundary_conditions_2D(settings, srh_all_Dict, nodeCoordinates)

#set up initial condition for for solution state variables for ODE solver
Q0 = hcat(h, q_x, q_y)

Q_ghost = hcat(h_ghostCells, q_x_ghostCells, q_y_ghostCells)

# Create the extra parameters struct
swe_extra_params = SWE2D_Extra_Parameters(
    active_param_name,
    settings,
    my_mesh_2D,
    nodeCoordinates,
    srh_all_Dict,
    boundary_conditions,
    swe_2D_constants,
    ManningN_cells,
    ManningN_ghostCells,
    inletQ_Length,
    inletQ_TotalQ,
    exitH_WSE,
    zb_cells,
    zb_ghostCells,
    zb_faces,
    S0
)

# Define the ODE function: p_extra is swe_params to pass more arguments to the ODE function
#@noinline 
# function swe_2d_ode(Q, p::Vector{T}, t::Float64, p_extra::SWEParameters{T}) where {T}
#     # p is just params_vector, use swe_params directly since it's in scope
#     dQdt = swe_2d_rhs(Q, p, t, p_extra)
#     return dQdt::Matrix{T}
# end

# Create the ODEFunction with the extra parameters struct passed in.
ode_f = ODEFunction((u, p, t) -> swe_2d_rhs(u, p, t, swe_extra_params); jac_prototype=jac_sparsity)

#########################
#Forward Simulation part#
#########################

if settings.bPerform_Forward_Simulation
    println("   Performing 2D SWE forward simulation ...")

    #perform forward simulation
    Hydrograd.swe_2D_forward_simulation(ode_f, Q0, params_vector, swe_extra_params, 
            zb_cells_truth, ManningN_values_truth, inlet_discharges_truth,
            case_path)
end


#########################
#Inversion part         # 
#########################

if settings.bPerform_Inversion

    println("   Performing inversion ...")

    #@show swe_2D_constants.tspan
    
    #perform inversion
    #@code_warntype 
    Hydrograd.swe_2D_inversion(ode_f, Q0, params_vector, swe_extra_params, case_path)

end


#########################   
#Sensitivity part       # 
#########################

if settings.bPerform_Sensitivity_Analysis

    println("   Performing sensitivity analysis ...")

    #perform sensitivity analysis
    Hydrograd.swe_2D_sensitivity(settings, my_mesh_2D, swe_2D_constants, ode_f, Q0, params_vector, active_range, param_ranges,
            nodeCoordinates, zb_cells, case_path)

end

#Timing 
end_time = now()  # Current date and time
elapsed_time = end_time - start_time
elapsed_seconds = Millisecond(elapsed_time).value / 1000
println("Elapsed time in seconds: $elapsed_seconds")

#restore the current directory
cd(current_dir)

println("All done!")
