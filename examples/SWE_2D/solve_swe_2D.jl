using Revise

using Dates

using JSON3

using JLD2

using IterTools

using PyCall

using Profile

using AdHydraulics

using OrdinaryDiffEq

using SparseArrays

#SciML
using SciMLSensitivity
using ComponentArrays

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

# Add some warnings/info for gradient computation
#@code_warntype(my_loss(ps))
SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true

#Timing 
start_time = now()  # Current date and time

include("process_SRH_2D_input.jl")
include("process_model_parameters_2D.jl")
include("process_ManningN_2D.jl")
include("process_BCs_2D.jl")
include("process_bed_2D.jl")
include("process_ICs_2D.jl")
include("semi_discretize_2D.jl")
include("misc_utilities_2D.jl")
include("custom_ODE_solver.jl")
include("process_inversion_results_2D.jl")

#print the banner
print_banner()

#directory to read from/save to (the same directory as the main file)
save_path = dirname(@__FILE__)

#read the control file
println("Reading control file...")
control_file = joinpath(save_path, "run_control.json")
control_dict = JSON3.read(open(control_file), Dict)

# Extract variables with their proper types
bVerbose = control_dict["bVerbose"]::Bool

# Control variables
srhhydro_file_name = control_dict["control_variables"]["srhhydro_file_name"]::String
bPerform_Forward_Simulation = control_dict["control_variables"]["bPerform_Forward_Simulation"]::Bool
bPerform_Inversion = control_dict["control_variables"]["bPerform_Inversion"]::Bool
bPerform_Sensitivity_Analysis = control_dict["control_variables"]["bPerform_Sensitivity_Analysis"]::Bool

# Currently only support one at a time: forward simulation, inversion, or sensitivity analysis
if bPerform_Forward_Simulation + bPerform_Inversion + bPerform_Sensitivity_Analysis != 1
    println("bPerform_Forward_Simulation = ", bPerform_Forward_Simulation)
    println("bPerform_Inversion = ", bPerform_Inversion)
    println("bPerform_Sensitivity_Analysis = ", bPerform_Sensitivity_Analysis)
    error("Currently only one at a time: forward simulation, inversion, or sensitivity analysis.")
end

# Time settings
bUse_srhhydro_time_settings = control_dict["time_settings"]["bUse_srhhydro_time_settings"]::Bool
tspan = Tuple(Float64.(control_dict["time_settings"]["tspan"]))  # Convert array to tuple
dt = Float64(control_dict["time_settings"]["dt"])
#dt_save = Float64(control_dict["time_settings"]["dt_save"])
#nSave = Int(control_dict["time_settings"]["nSave"])

# Forward simulation options
forward_simulation_solver = ""
forward_simulation_ode_solver = ""
forward_simulation_adaptive = false
forward_simulation_b_jac_sparsity = false
forward_simulation_nSave = 1
forward_simulation_initial_condition_options = ""
forward_simulation_initial_condition_file_name = ""
forward_simulation_initial_condition_values_from_file = nothing
forward_simulation_initial_condition_constant_values = []
forward_simulation_save_file_name = ""
forward_simulation_save_solution_truth_file_name = ""
    
if bPerform_Forward_Simulation
    forward_simulation_solver = control_dict["forward_simulation_options"]["forward_simulation_solver"]::String
    forward_simulation_ode_solver = control_dict["forward_simulation_options"]["forward_simulation_ode_solver"]::String
    forward_simulation_adaptive = control_dict["forward_simulation_options"]["forward_simulation_adaptive"]::Bool
    forward_simulation_b_jac_sparsity = control_dict["forward_simulation_options"]["forward_simulation_b_jac_sparsity"]::Bool
    forward_simulation_nSave = Int(control_dict["forward_simulation_options"]["forward_simulation_nSave"])
    forward_simulation_initial_condition_options = control_dict["forward_simulation_options"]["forward_simulation_initial_condition_options"]::String
    forward_simulation_initial_condition_file_name = control_dict["forward_simulation_options"]["forward_simulation_initial_condition_file_name"]::String
    forward_simulation_initial_condition_constant_values = Float64.(control_dict["forward_simulation_options"]["forward_simulation_initial_condition_constant_values"])

    if forward_simulation_initial_condition_options == "from_file"  #read the initial condition from a file
        forward_simulation_initial_condition_values_from_file = JSON3.read(open(joinpath(save_path, forward_simulation_initial_condition_file_name)), Dict)
    elseif forward_simulation_initial_condition_options == "constant"  #use constant initial condition
        forward_simulation_initial_condition_values_from_file = nothing
    else
        error("Invalid forward_simulation_initial_condition_options: $forward_simulation_initial_condition_options. Supported options: from_file, constant.")
    end

    forward_simulation_save_file_name = control_dict["forward_simulation_options"]["forward_simulation_save_file_name"]::String
    forward_simulation_save_solution_truth_file_name = control_dict["forward_simulation_options"]["forward_simulation_save_solution_truth_file_name"]::String
end

# Inversion options
active_param_names = []
inversion_parameter_initial_values_options = ""
inversion_parameter_initial_values_file_name = ""
inversion_parameter_initial_values_from_file = nothing
inversion_zb_initial_values = []
inversion_ManningN_initial_values = []
inversion_inlet_discharge_initial_values = []
inversion_bInversion_slope_loss = false
inversion_bInversion_u_loss = false
inversion_lower_bound_zb = 0.0
inversion_upper_bound_zb = 0.0
inversion_optimizer = ""
inversion_learning_rate = 0.0
inversion_max_iterations = 0
inversion_solution_truth_file_name = ""
inversion_forward_simulation_initial_condition_options = ""
inversion_forward_simulation_initial_condition_file_name = ""
inversion_forward_simulation_initial_condition_constant_values = []
inversion_forward_simulation_initial_condition_values_from_file = nothing
inversion_ode_solver = ""
inversion_ode_solver_adaptive = false
inversion_ode_solver_b_jac_sparsity = false

if bPerform_Inversion
    active_param_names = control_dict["inversion_options"]["active_param_names"]::Vector{String}
    inversion_parameter_initial_values_options = control_dict["inversion_options"]["inversion_parameter_initial_values_options"]::String
    inversion_parameter_initial_values_file_name = control_dict["inversion_options"]["inversion_parameter_initial_values_file_name"]::String
    inversion_zb_initial_values = control_dict["inversion_options"]["inversion_zb_initial_values"]::Vector{Float64}
    inversion_ManningN_initial_values = control_dict["inversion_options"]["inversion_ManningN_initial_values"]::Vector{Float64}
    inversion_inlet_discharge_initial_values = control_dict["inversion_options"]["inversion_inlet_discharge_initial_values"]::Vector{Float64}
    inversion_bInversion_slope_loss = control_dict["inversion_options"]["inversion_bInversion_slope_loss"]::Bool
    inversion_bInversion_u_loss = control_dict["inversion_options"]["inversion_bInversion_u_loss"]::Bool
    inversion_lower_bound_zb = control_dict["inversion_options"]["inversion_lower_bound_zb"]::Float64
    inversion_upper_bound_zb = control_dict["inversion_options"]["inversion_upper_bound_zb"]::Float64
    inversion_optimizer = control_dict["inversion_options"]["inversion_optimizer"]::String
    inversion_learning_rate = control_dict["inversion_options"]["inversion_learning_rate"]::Float64
    inversion_max_iterations = control_dict["inversion_options"]["inversion_max_iterations"]::Int
    inversion_solution_truth_file_name = control_dict["inversion_options"]["inversion_solution_truth_file_name"]::String
    inversion_forward_simulation_initial_condition_options = control_dict["inversion_options"]["inversion_forward_simulation_initial_condition_options"]::String
    inversion_forward_simulation_initial_condition_file_name = control_dict["inversion_options"]["inversion_forward_simulation_initial_condition_file_name"]::String
    inversion_forward_simulation_initial_condition_constant_values = Float64.(control_dict["inversion_options"]["inversion_forward_simulation_initial_condition_constant_values"])
    inversion_ode_solver = control_dict["inversion_options"]["ode_solver_options"]["ode_solver"]::String
    inversion_ode_solver_adaptive = control_dict["inversion_options"]["ode_solver_options"]["ode_solver_adaptive"]::Bool
    inversion_ode_solver_b_jac_sparsity = control_dict["inversion_options"]["ode_solver_options"]["ode_solver_b_jac_sparsity"]::Bool

    if inversion_parameter_initial_values_options == "from_file"
        inversion_parameter_initial_values_from_file = JSON3.read(open(joinpath(save_path, inversion_parameter_initial_values_file_name)), Dict)
    end
end

if bVerbose
    #print the control variables
    println("--------------------------------")

    println("Control variables:")
    println("   srhhydro_file_name = ", srhhydro_file_name)
    println("   bPerform_Forward_Simulation = ", bPerform_Forward_Simulation)
    println("   bPerform_Inversion = ", bPerform_Inversion)
    println("   bPerform_Sensitivity_Analysis = ", bPerform_Sensitivity_Analysis)

    println("Time settings (overwrite SRH-2D case if bUse_srhhydro_time_settings is false):")
    println("    bUse_srhhydro_time_settings = ", bUse_srhhydro_time_settings)
    println("    tspan = ", tspan)
    println("    dt = ", dt)

    if bPerform_Forward_Simulation
        println("Forward simulation options:")
        println("    solver = ", forward_simulation_solver)
        println("    ode_solver = ", forward_simulation_ode_solver)
        println("    adaptive = ", forward_simulation_adaptive)
        println("    jac_sparsity = ", forward_simulation_b_jac_sparsity)
        println("    nSave = ", forward_simulation_nSave)
        println("    initial_condition_options = ", forward_simulation_initial_condition_options)
        println("    initial_condition_file_name = ", forward_simulation_initial_condition_file_name)
        println("    save_file_name = ", forward_simulation_save_file_name)
        println("    save_solution_truth_file_name = ", forward_simulation_save_solution_truth_file_name)
    else
        println("No forward simulation is to be performed.")
    end

    if bPerform_Inversion
        println("Inversion options:")
        println("    active_params = ", active_param_names)
        println("    inversion_parameter_initial_values_options = ", inversion_parameter_initial_values_options)
        println("    inversion_parameter_initial_values_file_name = ", inversion_parameter_initial_values_file_name)

        if inversion_parameter_initial_values_options == "from_file"
            println("    inversion_parameter_initial_values_from_file = ", inversion_parameter_initial_values_from_file)
        elseif inversion_parameter_initial_values_options == "constant"
            println("    zb_initial_values = ", zb_initial_values)
            println("    ManningN_initial_values = ", ManningN_initial_values)
            println("    inlet_discharge_initial_values = ", inlet_discharge_initial_values)
        else
            error("Invalid inversion_parameter_initial_values_options: $inversion_parameter_initial_values_options. Supported options: from_file, constant.")
        end

        println("    slope_loss = ", bInversion_slope_loss)
        println("    u_loss = ", bInversion_u_loss)
        println("    zb_bounds = [", lower_bound_zb, ", ", upper_bound_zb, "]")
        println("    optimizer = ", optimizer)
        println("    learning_rate = ", learning_rate)
        println("    max_iterations = ", max_iterations)
        println("    inversion_solution_truth_file_name = ", inversion_solution_truth_file_name)
        println("    ode_solver = ", inversion_ode_solver)
        println("    ode_solver_adaptive = ", inversion_ode_solver_adaptive)
        println("    ode_solver_b_jac_sparsity = ", inversion_ode_solver_b_jac_sparsity)
    else
        println("No inversion is to be performed.")
    end

    if bPerform_Sensitivity_Analysis
        println("Sensitivity analysis is to be performed.")
    else
        println("No sensitivity analysis is to be performed.")
    end

    println("--------------------------------")
end

println("Finished reading control file.\n")

#define a swe_2D_constants object with some values from the control file
swe_2D_constants = swe_2D_consts(t=tspan[1], dt=dt, tStart=tspan[1], tEnd=tspan[2])

#read data from SRH-2D hydro, geom, and material files; it aslo create the 2D mesh.
srh_all_Dict = process_SRH_2D_input(srhhydro_file_name)

#update swe_2D_constants based on the SRH-2D data
if bUse_srhhydro_time_settings
    update_swe_2D_constants!(swe_2D_constants, srh_all_Dict)
end

#Get the 2D mesh 
my_mesh_2D = srh_all_Dict["my_mesh_2D"]

#mesh related variables which might be mutable (should not be in struct mesh_2D)
#nodeCoordinates: Float64 2D array [numOfNodes, 3]
nodeCoordinates = srh_all_Dict["srhgeom_obj"].nodeCoordinates

#  setup bed elevation 
#If performing inversion, set the bed elevation to zero to start with
if bPerform_Inversion
    nodeCoordinates[:,3] .= 0.0
end

#setup bed elevation: computer zb at cell centers (zb_cells) from nodes, 
#then interpolate zb from cell to face (zb_faces) and ghost cells (zb_ghostCells),
#and compute the bed slope at cells (S0)
zb_cells, zb_ghostCells, zb_faces, S0 = setup_bed(my_mesh_2D, nodeCoordinates, true)

zb_cells_truth = zeros(size(zb_cells))

#If performing forward simulation, make a copy of the bathymetry truth: zb_cell_truth; otherwise, it is zero.
if bPerform_Forward_Simulation
    zb_cells_truth = deepcopy(zb_cells)
end

#get the true Manning's n and inlet discharges
srhhydro_ManningsN_Dict = srh_all_Dict["srhhydro_ManningsN"]

ManningN_values_truth = Float64[srhhydro_ManningsN_Dict[i] for i in 0:(length(srhhydro_ManningsN_Dict)-1)]

#get the true inlet discharges
srhhydro_inletQ_Dict = srh_all_Dict["srhhydro_IQParams"]

inlet_discharges_truth = [parse(Float64, srhhydro_inletQ_Dict[i][1]) for i in 1:(length(srhhydro_inletQ_Dict))]

if bPerform_Forward_Simulation && bVerbose
    println("True zb_cells = ", zb_cells_truth)
    println("True Manning's n values = ", ManningN_values_truth)
    println("True inlet discharges = ", inlet_discharges_truth)
end

#define the parameters array (nothing for forward simulation)
zb_cells_param = nothing
ManningN_list_param = nothing
inlet_discharges_param = nothing

if bPerform_Inversion || bPerform_Sensitivity_Analysis
    if inversion_parameter_initial_values_options == "constant"
        #bed elevation
        zb_cells_param = zeros(zb_initial_values[1], my_mesh_2D.numOfCells)  #only use the first value of zb_initial_values even it is a vector.

        #Manning's n
        ManningN_list_param = deepcopy(ManningN_initial_values)   
        
        #inlet discharges
        inlet_discharges_param = deepcopy(inlet_discharges_initial_values)  
    elseif inversion_parameter_initial_values_options == "from_file"
        zb_cells_param = inversion_parameter_initial_values_from_file["zb_cells_param"]
        ManningN_list_param = inversion_parameter_initial_values_from_file["ManningN_list_param"]
        inlet_discharges_param = inversion_parameter_initial_values_from_file["inlet_discharges_param"]
    else
        error("Invalid zb_initial_values_options: $zb_initial_values_options. Supported options: zero, from_file.")
    end

    #consistency check
    #make sure the length of zb_cells_param is the same as the number of cells
    if length(zb_cells_param) != my_mesh_2D.numOfCells
        error("The length of zb_cells_param is not the same as the number of cells.")
    end

    #make sure the length of ManningN_list_param is the same as the number of materials
    srhmat_numOfMaterials = srh_all_Dict["srhmat_numOfMaterials"]
    if length(ManningN_list_param) != srhmat_numOfMaterials
        error("The length of ManningN_list_param is not the same as the number of materials.")
    end

    #make sure the length of inlet_discharges_param is the same as the number of inletQ_BCs
    nInletQ_BCs = srh_all_Dict["nInletQ_BCs"]
    if length(inlet_discharges_param) != nInletQ_BCs
        error("The length of inlet_discharges_param is not the same as the number of inletQ_BCs.")
    end
end

#preprocess: create a ComponentArray for the model parameters for 2D shallow water equations
params_array, active_params = preprocess_model_parameters_2D(bPerform_Forward_Simulation, bPerform_Inversion, bPerform_Sensitivity_Analysis, zb_cells_param, ManningN_list_param, inlet_discharges_param, active_param_names)

#Initial setup of Manning's n using the SRH-2D data (if performing inversion on Manning's n, ManningN_cells will be updated later in the inversion process)
ManningN_cells, ManningN_ghostCells = setup_ManningN(my_mesh_2D, srh_all_Dict)

# setup initial conditions 
eta = zeros(my_mesh_2D.numOfCells)          #free surface elevation at cells 
h = zeros(my_mesh_2D.numOfCells)            #water depth at cells 
q_x = zeros(my_mesh_2D.numOfCells)          #q_x=hu at cells 
q_y = zeros(my_mesh_2D.numOfCells)          #q_y=hv at cells 

eta_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #free surface elevation at ghost cells 
h_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)            #water depth at ghost cells 
q_x_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #q_x=hu at ghost cells 
q_y_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #q_y=hv at ghost cells 

total_water_volume = []   #total volume of water in the domain 

#setup initial condition for eta, h, q_x, q_y
initial_condition_options = nothing
initial_condition_constant_values = nothing
initial_condition_values_from_file = nothing

if bPerform_Forward_Simulation
    initial_condition_options = forward_simulation_initial_condition_options
    initial_condition_constant_values = forward_simulation_initial_condition_constant_values
    initial_condition_values_from_file = forward_simulation_initial_condition_values_from_file
elseif bPerform_Inversion 
    initial_condition_options = inversion_forward_simulation_initial_condition_options
    initial_condition_constant_values = inversion_forward_simulation_initial_condition_constant_values
    initial_condition_values_from_file = inversion_forward_simulation_initial_condition_values_from_file
elseif bPerform_Sensitivity_Analysis
    #not implemented yet
    error("Sensitivity analysis is not implemented yet.")
else
    error("Invalid bPerform_Forward_Simulation, bPerform_Inversion, bPerform_Sensitivity_Analysis. No initial condition is to be setup.")
end 

setup_initial_condition!(initial_condition_options, initial_condition_constant_values, initial_condition_values_from_file, my_mesh_2D, nodeCoordinates, eta, zb_cells, h, q_x, q_y, swe_2D_constants, true)

#setup ghost cells for initial condition
setup_ghost_cells_initial_condition!(my_mesh_2D, eta, h, q_x, q_y, eta_ghostCells, h_ghostCells, q_x_ghostCells, q_y_ghostCells)

#create and preprocess boundary conditions: boundary_conditions only contains the static information of the boundaries.
boundary_conditions, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length, inletQ_TotalA, inletQ_DryWet, 
exitH_WSE, exitH_H, exitH_A, wall_H, wall_A, symm_H, symm_A = initialize_boundary_conditions_2D(srh_all_Dict, nodeCoordinates)

#set up initial condition for ODE solver
Q0 = hcat(h, q_x, q_y)  

Q_ghost = hcat(h_ghostCells, q_x_ghostCells, q_y_ghostCells)   #ghost cell values of Q

# Define the ODE function 
function swe_2d_ode(Q, params_array, t)

    dQdt = swe_2d_rhs(Q, params_array, t, bPerform_Forward_Simulation, bPerform_Inversion, bPerform_Sensitivity_Analysis, 
                      my_mesh_2D, boundary_conditions, swe_2D_constants, ManningN_cells, ManningN_ghostCells, inletQ_Length, inletQ_TotalQ, exitH_WSE,
                      zb_cells, zb_ghostCells, zb_faces, S0,
                      active_param_names)    

    return dQdt
end

# Create the ODEFunction with the typed function
ode_f = nothing

jac_sparsity = nothing

if (bPerform_Forward_Simulation && forward_simulation_b_jac_sparsity) || 
    (bPerform_Inversion && inversion_ode_solver_b_jac_sparsity) ||
    (bPerform_Sensitivity_Analysis && sensitivity_analysis_ode_solver_b_jac_sparsity)

    # Populate the Jacobian sparsity pattern
    # Assume each variable depends on itself and its neighbors
    # we have 3 variables (h, q_x, q_y) for each cell
    jac_sparsity = spzeros(3*my_mesh_2D.numOfCells, 3*my_mesh_2D.numOfCells)

    for cellID in 1:my_mesh_2D.numOfCells
        # Self-dependence (diagonal entries)
        jac_sparsity[3*(cellID-1)+1, 3*(cellID-1)+1] = 1.0  # h -> h
        jac_sparsity[3*(cellID-1)+2, 3*(cellID-1)+2] = 1.0  # hu -> hu
        jac_sparsity[3*(cellID-1)+3, 3*(cellID-1)+3] = 1.0  # hv -> hv
        
        # Neighbor-dependence
        for neighbor in my_mesh_2D.cellNeighbors_Dict[cellID]

            #only need interior neighbors
            if neighbor > 0
                jac_sparsity[3*(cellID-1)+1, 3*(neighbor-1)+1] = 1.0  # h -> h (neighbor)
                jac_sparsity[3*(cellID-1)+2, 3*(neighbor-1)+2] = 1.0  # hu -> hu (neighbor)
                jac_sparsity[3*(cellID-1)+3, 3*(neighbor-1)+3] = 1.0  # hv -> hv (neighbor)
            end
        end
    end
end

ode_f = ODEFunction(swe_2d_ode; jac_prototype=jac_sparsity)


#########################
#Forward Simulation part#
#########################

if bPerform_Forward_Simulation
    println("   Performing 2D SWE forward simulation ...")

    # time information 
    tspan = (swe_2D_constants.tStart, swe_2D_constants.tEnd)
    dt = swe_2D_constants.dt
    t = tspan[1]:dt:tspan[2]

    #define the time for saving the results
    dt_save = (tspan[2] - tspan[1])/forward_simulation_nSave
    t_save = tspan[1]:dt_save:tspan[2]
    println("t_save = ", t_save)

    #define the ODE problem
    prob = ODEProblem(ode_f, Q0, tspan, params_array)
 
    if forward_simulation_solver == "SciML"

        println("       with SciML solver ...")
    
        if forward_simulation_ode_solver == "Tsit5()"
            sol = solve(prob, Tsit5(), adaptive=forward_simulation_adaptive, dt=dt, saveat=t_save)
        else
            sol = solve(prob, Tsit5(), adaptive=true, dt=dt, saveat=t_save)
        end

        # #save the simulation solution results
        jldsave(joinpath(save_path, forward_simulation_save_file_name); sol)

        #save the simulation results (h, u, v) at the last time step to a json file (to be used as ground truth for inversion)
        h_truth = Array(sol)[:,1,end]
        eta_truth = h_truth .+ zb_cells_truth
        u_truth = Array(sol)[:,2,end]./Array(sol)[:,1,end]
        v_truth = Array(sol)[:,3,end]./Array(sol)[:,1,end]
        open(joinpath(save_path, forward_simulation_save_solution_truth_file_name), "w") do io
            JSON3.pretty(io, Dict("eta_truth" => eta_truth, "h_truth" => h_truth, "u_truth" => u_truth, "v_truth" => v_truth, "zb_cells_truth" => zb_cells_truth, "ManningN_values_truth" => ManningN_values_truth, "inlet_discharges_truth" => inlet_discharges_truth))
            println(io)
        end

        swe_2D_save_results_SciML(sol, total_water_volume, my_mesh_2D, nodeCoordinates, zb_cells, save_path)
    
    elseif solver_choice == "customized"    #My own ODE Solver

        println("   Performing 2D SWE simulation with MyOwn solver ...")
        
        sol = my_solve(para, Q0, my_mesh_2D, tspan, dt)

        # #save the results
        # #save the simulation solution results
        jldsave(joinpath(save_path, "simulation_solution.jld2"); sol)

        swe_2D_save_results_custom(sol, total_water_volume, my_mesh_2D, zb_cells, save_path)
    else
        println("   Wrong solver choice. Supported solvers: SciML, customized. No forward simulation is performed.")
    end
end


#########################
#Inversion part         # 
#########################

if bPerform_Inversion

    println("   Performing inversion ...")

    #open the forward simulation result (as the ground truth)
    sol_truth = JSON3.read(open(joinpath(save_path, inversion_solution_truth_file_name)), Dict)

    WSE_truth = sol_truth["eta_truth"]
    h_truth = sol_truth["h_truth"]
    u_truth = sol_truth["u_truth"]
    v_truth = sol_truth["v_truth"]
    zb_cells_truth = sol_truth["zb_cells_truth"]
    ManningN_values_truth = sol_truth["ManningN_values_truth"]
    inlet_discharges_truth = sol_truth["inlet_discharges_truth"]

    #inversion parameters: 
    # initial guess for the inversion parameters
   
    params_array_init = params_array
    
    
    function predict(θ)
        #Zygote.ignore() do
        #    println("Current parameters: ", ForwardDiff.value.(θ))
        #end

        #See https://docs.sciml.ai/SciMLSensitivity/dev/faq/ for the choice of AD type (sensealg)
        #If not specified, the default is a smart polyalgorithm used to automatically determine the most appropriate method for a given equation.
        #Some options are:
        # 1. BacksolveAdjoint(autojacvec=ZygoteVJP())
        # 2. InterpolatingAdjoint(autojacvec=ZygoteVJP())
        # 3. ReverseDiffAdjoint(autojacvec=ZygoteVJP())
        # 4. TrackerAdjoint(autojacvec=ZygoteVJP())
        # 5. ZygoteVJP()
        # 6. TrackerVJP()
        # 7. ReverseDiffVJP()
        # 8. InterpolatingVJP()
        # 9. BacksolveVJP()
        #Array(solve(prob, Heun(), adaptive=false, p=θ, dt=dt, saveat=t))[:,1,end]
        #sol = solve(prob, Tsit5(), adaptive=false, p=θ, dt=dt, saveat=t_save)  #[:,1,end]  #works, but takes a long time (10x slower than with sensealg=ForwardDiffSensitivity()). Due to adaptive=false?
        sol = solve(prob, Tsit5(), adaptive=true, p=θ, dt=dt, saveat=t_save) #adaptive=true is default. This is 10x faster than with adaptive=false.
        #sol = solve(prob, Euler(), adaptive=false, p=θ, dt=dt, saveat=t_save)  #[:,1,end]

        #sol = solve(prob, Tsit5(), p=θ, dt=dt, saveat=t_save; sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP())) #only works for short time span
        #sol = solve(prob, Tsit5(), p=θ, dt=dt, saveat=t_save; sensealg=ForwardDiffSensitivity())   #runs, but very slow
        #sol = solve(prob, Tsit5(), adaptive=false, p=θ, dt=dt, saveat=t_save; sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()))
        #sol = solve(prob, Tsit5(), adaptive=false, p=θ, dt=dt, saveat=t_save; sensealg=GaussAdjoint(autojacvec=ZygoteVJP()))  #not working

        #if !SciMLBase.successful_retcode(sol)
        #    # Return a high cost instead of NaN
        #    return fill(convert(eltype(θ), Inf), size(sol[1]))
        #end
        
        return sol
    end


    ## Define Loss function
    function loss(θ)
        pred = predict(θ)            #Forward prediction with current parameter θ 
                
        #\theta is the current bed elevation at cells
        l = pred[:,1,end] .+ θ .- WSE_truth  #loss = free surface elevation mismatch

        loss_pred_eta = sum(abs2, l)

        loss_pred_uv = zero(eltype(θ))

        #  # Add small epsilon to prevent division by zero
        # ϵ = sqrt(eps(eltype(θ)))
        
        # if bInversion_include_u      #if also include u in the loss 
        #     l_u = pred[:,2,end]./(pred[:,1,end] .+ ϵ) .- u_truth
        #     l_v = pred[:,3,end]./(pred[:,1,end] .+ ϵ) .- v_truth

        #     loss_pred_uv = sum(abs2, l_u) + sum(abs2, l_v)
        # end 

        loss_pred = loss_pred_eta + loss_pred_uv

        loss_slope = zero(eltype(θ))

        # if bInversion_slope_loss    #if bed slope is included in the loss 
        #     loss_slope = calc_slope_loss(θ, my_mesh_2D)
        # end 

        loss_total = loss_pred + loss_slope

        return loss_total, loss_pred, loss_pred_eta, loss_pred_uv, loss_slope, pred
        # return loss_total

        #Zygote.ignore() do
        #    println("loss_pred_eta = ", loss_pred_eta)
        #end

        #return loss_pred_eta
    end


    LOSS = []                              # Loss accumulator
    PRED = []                              # prediction accumulator
    PARS = []                              # parameters accumulator

    callback = function (θ, loss_total, loss_pred, loss_pred_eta, loss_pred_uv, loss_slope, pred) #callback function to observe training
        iter = size(LOSS)[1]  #get the inversion iteration number (=length of LOSS array)
        println("      iter, loss_total, loss_pred, loss_pred_eta, loss_pred_uv, loss_slope = ", iter, ", ", 
                  loss_total, ", ", loss_pred, ", ", loss_pred_eta, ", ", loss_pred_uv, ", ", loss_slope)

        append!(PRED, [pred[:,1,end]])
        append!(LOSS, [[loss_total, loss_pred, loss_pred_eta, loss_pred_uv, loss_slope, pred]])

        if !isa(θ, Vector{Float64})  #NLopt returns an optimization object, not an arrary
            #println("theta.u = ", θ.u)
            append!(PARS, [copy(θ.u)])
        else
            append!(PARS, [θ])
        end

        #if l > 1e-9
        #    false
        #else 
        #    true   #force the optimizer to stop 
        #end

        false
    end

    # Define AD type choice for optimization's gradient computation
    #The following is from SciMLSensitivity documentation regarding the choice of AD 
    #  AutoForwardDiff(): The fastest choice for small optimizations
    #  AutoReverseDiff(compile=false): A fast choice for large scalar optimizations
    #  AutoTracker(): Like ReverseDiff but GPU-compatible
    #  AutoZygote(): The fastest choice for non-mutating array-based (BLAS) functions
    #  AutoFiniteDiff(): Finite differencing, not optimal but always applicable
    #  AutoModelingToolkit(): The fastest choice for large scalar optimizations
    #  AutoEnzyme(): Highly performant AD choice for type stable and optimized code

    adtype = Optimization.AutoZygote()         #
    #adtype = Optimization.AutoForwardDiff()   #ForwardDiff.jl is not supported for this problem.
    #adtype = Optimization.AutoReverseDiff()
    #adtype = Optimization.AutoEnzyme()         #

    # Define the optimization function
    #From SciMLSensitivity documentation: https://docs.sciml.ai/Optimization/stable/API/optimization_function/
    # OptimizationFunction{iip}(f, adtype::AbstractADType = NoAD();
    #                       grad = nothing, hess = nothing, hv = nothing,
    #                       cons = nothing, cons_j = nothing, cons_jvp = nothing,
    #                       cons_vjp = nothing, cons_h = nothing,
    #                       hess_prototype = nothing,
    #                       cons_jac_prototype = nothing,
    #                       cons_hess_prototype = nothing,
    #                       observed = __has_observed(f) ? f.observed : DEFAULT_OBSERVED_NO_TIME,
    #                       lag_h = nothing,
    #                       hess_colorvec = __has_colorvec(f) ? f.colorvec : nothing,
    #                       cons_jac_colorvec = __has_colorvec(f) ? f.colorvec : nothing,
    #                       cons_hess_colorvec = __has_colorvec(f) ? f.colorvec : nothing,
    #                       lag_hess_colorvec = nothing,
    #                       sys = __has_sys(f) ? f.sys : nothing)

    optf = Optimization.OptimizationFunction((θ, p) -> loss(θ), adtype)   #\theta is the parameter to be optimized; p is not used.

    # Define the bounds for the parameter (only applicable for some optimizers which support lb and ub)
    lb_p = zeros(my_mesh_2D.numOfCells)
    lb_p .= -0.1  

    ub_p = zeros(my_mesh_2D.numOfCells)
    ub_p .= 0.3 

    # Define the optimization problem
    #From SciMLSensitivity documentation: https://docs.sciml.ai/Optimization/stable/API/optimization_problem/
    #OptimizationProblem{iip}(f, u0, p = SciMLBase.NullParameters(),;
    #                          lb = nothing,
    #                          ub = nothing,
    #                          lcons = nothing,
    #                          ucons = nothing,
    #                          sense = nothing,
    #                          kwargs...)

    #optprob = Optimization.OptimizationProblem(optf, ps, lb=lb_p, ub=ub_p)
    optprob = Optimization.OptimizationProblem(optf, ps)

    # Solve the optimization problem
    #From SciMLSensitivity documentation: https://docs.sciml.ai/Optimization/stable/API/optimization_solution/
    #Returned optimization solution Fields:
    # u: the representation of the optimization's solution.
    # cache::AbstractOptimizationCache: the optimization cache` that was solved.
    # alg: the algorithm type used by the solver.
    # objective: Objective value of the solution
    # retcode: the return code from the solver. Used to determine whether the solver solved successfully or whether it exited due to an error. For more details, see the return code documentation.
    # original: if the solver is wrapped from a external solver, e.g. Optim.jl, then this is the original return from said solver library.
    # stats: statistics of the solver, such as the number of function evaluations required.
    
    #res = Optimization.solve(optprob, PolyOpt(), callback = callback)  #PolyOpt does not support lb and ub 
    #res = Optimization.solve(optprob, NLopt.LD_LBFGS(), callback = callback)   #very fast 
    #res = Optimization.solve(optprob, Optim.BFGS(), callback=callback; iterations=30, maxiters=40, f_calls_limit=20, show_trace=true)
    #res = Optimization.solve(optprob, Optim.BFGS(), callback=callback, maxiters = 100; show_trace=false)  #f_tol=1e-3, iterations=10, local_maxiters = 10
    #res = Optimization.solve(optprob, Optim.LBFGS(), callback=callback)  #oscilates around 1e-7
    #res = Optimization.solve(optprob, Optim.Newton(), callback=callback)  #error: not supported as the Fminbox optimizer
    #res = Optimization.solve(optprob, Optim.GradientDescent(), callback=callback)  #very slow decrease in loss 
    res = Optimization.solve(optprob, Adam(0.01), callback=callback, maxiters=100)

    #@time res = Optimization.solve(optprob, Adam(0.01), callback=callback, maxiters=100)
    #@profile res = Optimization.solve(optprob, Adam(0.01), callback=callback, maxiters=100)
    #Profile.print()
    
    #@show res
    

    #save the inversion results
    jldsave(joinpath(save_path, "inversion_results.jld2"); LOSS, PRED, PARS)

    #process the inversion results
    println("   Processing inversion results ...")

    #open the simulation result (as the ground truth)
    sol_truth = load(joinpath(save_path, "simulation_solution.jld2"))["sol"]
    h_truth = Array(sol_truth)[:,1,end]
    u_truth = Array(sol_truth)[:,2,end]./Array(sol_truth)[:,1,end]
    v_truth = Array(sol_truth)[:,3,end]./Array(sol_truth)[:,1,end]

    zb_cells_truth = load(joinpath(save_path, "zb_cells_truth.jld2"))["zb_cells_truth"]

    WSE_truth = h_truth .+ zb_cells_truth

    #process inversion results
    inversion_results_file_name = joinpath(save_path, "inversion_results.jld2")
    process_inversion_results(inversion_results_file_name, my_mesh_2D, nodeCoordinates, zb_cells_truth, h_truth, u_truth, v_truth, WSE_truth)
end 


#Timing 
end_time = now()  # Current date and time
elapsed_time = end_time - start_time
elapsed_seconds = Millisecond(elapsed_time).value / 1000
println("Elapsed time in seconds: $elapsed_seconds")

println("All done!")
