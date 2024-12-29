using Revise

using Dates

using JSON3

using JLD2

using IterTools

using PyCall

using Profile

using Hydrograd

using SparseArrays

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

# Add some warnings/info for gradient computation
SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true

#Timing 
start_time = now()  # Current date and time

#include the files to define the variables
include("define_variables_2D.jl")

#print the banner
print_banner()

#directory to read from/save to (the same directory as this main file)
save_path = dirname(@__FILE__)

# Read and parse control file
println("Reading control file...")
control_file = joinpath(save_path, "run_control.json")
settings = parse_control_file(control_file)

#define a swe_2D_constants object with some values from the control file
swe_2D_constants = swe_2D_consts(t=settings.time_settings.tspan[1], dt=settings.time_settings.dt, tStart=settings.time_settings.tspan[1], tEnd=settings.time_settings.tspan[2])

#read data from SRH-2D hydro, geom, and material files; it aslo create the 2D mesh.
srh_all_Dict = process_SRH_2D_input(settings)

#update swe_2D_constants based on the SRH-2D data
if settings.time_settings.bUse_srhhydro_time_settings
    update_swe_2D_constants!(swe_2D_constants, srh_all_Dict)
end

#Get the 2D mesh 
my_mesh_2D = srh_all_Dict["my_mesh_2D"]

#mesh related variables which might be mutable (should not be in struct mesh_2D)
#nodeCoordinates: Float64 2D array [numOfNodes, 3]
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
zb_cells, zb_ghostCells, zb_faces, S0 = setup_bed(settings, my_mesh_2D, nodeCoordinates, true)

zb_cells_truth = zeros(size(zb_cells))

#If performing forward simulation, make a copy of the bathymetry truth: zb_cell_truth; otherwise, it is zero.
if settings.bPerform_Forward_Simulation
    zb_cells_truth = deepcopy(zb_cells)
end

#get the true Manning's n and inlet discharges
srhhydro_ManningsN_Dict = srh_all_Dict["srhhydro_ManningsN"]

ManningN_values_truth = Float64[srhhydro_ManningsN_Dict[i] for i in 0:(length(srhhydro_ManningsN_Dict)-1)]

#get the true inlet discharges (could be nothing if no inletQ_BCs)
srhhydro_inletQ_Dict = srh_all_Dict["srhhydro_IQParams"]

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

#define the parameters array (nothing for forward simulation)
zb_cells_param = nothing
ManningN_list_param = nothing
inlet_discharges_param = nothing

#for forward simulation, the parameter values are from SRH-2D data
if settings.bPerform_Forward_Simulation
    zb_cells_param = deepcopy(zb_cells_truth)
    ManningN_list_param = deepcopy(ManningN_values_truth)

    if !isnothing(inlet_discharges_truth)
        inlet_discharges_param = deepcopy(inlet_discharges_truth)
    end
end

if settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis
    if settings.inversion_settings.parameter_initial_values_options == "constant"
        #bed elevation
        zb_cells_param = fill(settings.inversion_settings.zb_initial_values[1], my_mesh_2D.numOfCells)

        #Manning's n
        ManningN_list_param = deepcopy(settings.inversion_settings.ManningN_initial_values)

        #inlet discharges
        if !isnothing(inlet_discharges_truth)
            inlet_discharges_param = deepcopy(settings.inversion_settings.inlet_discharge_initial_values)
        end
    elseif settings.inversion_settings.parameter_initial_values_options == "from_file"
        zb_cells_param = settings.inversion_settings.parameter_initial_values_from_file["zb_cells_param"]
        ManningN_list_param = settings.inversion_settings.parameter_initial_values_from_file["ManningN_list_param"]

        if !isnothing(inlet_discharges_truth)
            inlet_discharges_param = settings.inversion_settings.parameter_initial_values_from_file["inlet_discharges_param"]
        end
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

    if !isnothing(inlet_discharges_truth)
        if length(inlet_discharges_param) != nInletQ_BCs
            error("Length mismatch: inlet_discharges_param ($(length(inlet_discharges_param))) != nInletQ_BCs ($nInletQ_BCs)")
        end
    else    #if inlet_discharges_truth is nothing, then nInletQ_BCs should be 0
        if nInletQ_BCs > 0
            error("No inletQ_BCs are defined in the SRH-2D data, but nInletQ_BCs is greater than 0. Please check the SRH-2D data.")
        end
    end
end

#preprocess: create a array of model parameters for 2D shallow water equations
#params_array: the 1D array of all parameters (zb_cells_param, ManningN_list_param, inlet_discharges_param)
#active_range: the range of active parameters 
#param_ranges: the range of each parameter
params_array, active_range, param_ranges = preprocess_model_parameters_2D(settings, zb_cells_param, ManningN_list_param, inlet_discharges_param)

#Initial setup of Manning's n using the SRH-2D data (if performing inversion on Manning's n, ManningN_cells will be updated later in the inversion process)
ManningN_cells, ManningN_ghostCells = setup_ManningN(settings, my_mesh_2D, srh_all_Dict)

# setup initial conditions 
wse = zeros(my_mesh_2D.numOfCells)          #free surface elevation at cells 
h = zeros(my_mesh_2D.numOfCells)            #water depth at cells 
q_x = zeros(my_mesh_2D.numOfCells)          #q_x=hu at cells 
q_y = zeros(my_mesh_2D.numOfCells)          #q_y=hv at cells 

wse_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #free surface elevation at ghost cells 
h_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)            #water depth at ghost cells 
q_x_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #q_x=hu at ghost cells 
q_y_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #q_y=hv at ghost cells 

total_water_volume = []   #total volume of water in the domain 

#setup initial condition for wse, h, q_x, q_y
initial_condition_options = nothing
initial_condition_constant_values = nothing
initial_condition_values_from_file = nothing

if settings.bPerform_Forward_Simulation
    initial_condition_options = settings.forward_settings.initial_condition_options
    initial_condition_constant_values = settings.forward_settings.initial_condition_constant_values
    initial_condition_values_from_file = settings.forward_settings.initial_condition_values_from_file
elseif settings.bPerform_Inversion
    initial_condition_options = settings.inversion_settings.forward_simulation_initial_condition_options
    initial_condition_constant_values = settings.inversion_settings.forward_simulation_initial_condition_constant_values
    initial_condition_values_from_file = settings.inversion_settings.forward_simulation_initial_condition_values_from_file
elseif settings.bPerform_Sensitivity_Analysis
    #not implemented yet
    error("Sensitivity analysis is not implemented yet.")
else
    error("Invalid bPerform_Forward_Simulation, bPerform_Inversion, bPerform_Sensitivity_Analysis. No initial condition is to be setup.")
end

setup_initial_condition!(settings, initial_condition_options, initial_condition_constant_values, initial_condition_values_from_file,
    my_mesh_2D, nodeCoordinates, wse, zb_cells, h, q_x, q_y, swe_2D_constants, true)

#setup ghost cells for initial condition
setup_ghost_cells_initial_condition!(settings, my_mesh_2D, wse, h, q_x, q_y, wse_ghostCells, h_ghostCells, q_x_ghostCells, q_y_ghostCells)

#create and preprocess boundary conditions: boundary_conditions only contains the static information of the boundaries.
boundary_conditions, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length, inletQ_TotalA, inletQ_DryWet,
exitH_WSE, exitH_H, exitH_A, wall_H, wall_A, symm_H, symm_A = initialize_boundary_conditions_2D(settings, srh_all_Dict, nodeCoordinates)

#set up initial condition for ODE solver
Q0 = hcat(h, q_x, q_y)

Q_ghost = hcat(h_ghostCells, q_x_ghostCells, q_y_ghostCells)   #ghost cell values of Q

# Define the ODE function 
function swe_2d_ode(Q, params_array, t)

    dQdt = swe_2d_rhs(Q, params_array, active_range, param_ranges, t, settings,
        my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants, ManningN_cells, ManningN_ghostCells, inletQ_Length, inletQ_TotalQ, exitH_WSE,
        zb_cells, zb_ghostCells, zb_faces, S0)

    return dQdt
end

# Create the ODEFunction with the typed function
ode_f = nothing

jac_sparsity = nothing

if (settings.bPerform_Forward_Simulation && settings.forward_settings.ode_solver_b_jac_sparsity) ||
   (settings.bPerform_Inversion && settings.inversion_settings.ode_solver_b_jac_sparsity) ||
   (settings.bPerform_Sensitivity_Analysis && settings.sensitivity_analysis_settings.ode_solver_b_jac_sparsity)

    # Populate the Jacobian sparsity pattern
    # Assume each variable depends on itself and its neighbors
    # we have 3 variables (h, q_x, q_y) for each cell
    jac_sparsity = spzeros(3 * my_mesh_2D.numOfCells, 3 * my_mesh_2D.numOfCells)

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

#define the ODE function
ode_f = ODEFunction(swe_2d_ode; jac_prototype=jac_sparsity)

# time information (the same for forward simulation, inversion, and sensitivity analysis)
tspan = (swe_2D_constants.tStart, swe_2D_constants.tEnd)
dt = swe_2D_constants.dt
t = tspan[1]:dt:tspan[2]

#########################
#Forward Simulation part#
#########################

if settings.bPerform_Forward_Simulation
    println("   Performing 2D SWE forward simulation ...")

    #define the time for saving the results (for forward simulation, which may be different from the time for saving the results in inversion)
    dt_save = (tspan[2] - tspan[1]) / settings.forward_settings.nSave
    t_save = tspan[1]:dt_save:tspan[2]

    if settings.bVerbose
        println("t_save = ", t_save)
    end

    #define the ODE problem
    prob = ODEProblem(ode_f, Q0, tspan, params_array)

    if settings.forward_settings.solver == "SciML"

        println("       with SciML solver ...")

        if settings.forward_settings.ode_solver == "Tsit5()"
            sol = solve(prob, Tsit5(), adaptive=settings.forward_settings.ode_solver_adaptive, dt=dt, saveat=t_save)
        else
            sol = solve(prob, Tsit5(), adaptive=true, dt=dt, saveat=t_save)
        end

        # #save the simulation solution results
        jldsave(joinpath(save_path, settings.forward_settings.save_file_name); sol)

        #save the simulation results (h, u, v) at the last time step to a json file (to be used as ground truth for inversion)
        h_truth = Array(sol)[:, 1, end]
        wse_truth = h_truth .+ zb_cells_truth
        u_truth = Array(sol)[:, 2, end] ./ Array(sol)[:, 1, end]
        v_truth = Array(sol)[:, 3, end] ./ Array(sol)[:, 1, end]
        open(joinpath(save_path, settings.forward_settings.save_solution_truth_file_name), "w") do io
            JSON3.pretty(io, Dict("wse_truth" => wse_truth, "h_truth" => h_truth, "u_truth" => u_truth, "v_truth" => v_truth, "zb_cells_truth" => zb_cells_truth, "ManningN_values_truth" => ManningN_values_truth, "inlet_discharges_truth" => inlet_discharges_truth))
            println(io)
        end

        swe_2D_save_results_SciML(sol, total_water_volume, my_mesh_2D, nodeCoordinates, zb_cells, save_path)

    elseif settings.forward_settings.solver == "customized"    #My own ODE Solver

        println("   Performing 2D SWE simulation with MyOwn solver ...")

        sol = custom_ODE_solve(params_array, Q0, my_mesh_2D, tspan, dt)

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

if settings.bPerform_Inversion

    println("   Performing inversion ...")

    #define the time for saving the results (for inversion, which may be different from the time for saving the results in forward simulation)
    dt_save = (tspan[2] - tspan[1]) / settings.inversion_settings.ode_solver_nSave
    t_save = tspan[1]:dt_save:tspan[2]

    if settings.bVerbose
        println("t_save = ", t_save)
    end

    #open the forward simulation result (as the ground truth)
    sol_truth = JSON3.read(open(joinpath(save_path, settings.inversion_settings.inversion_truth_file_name)), Dict)

    WSE_truth = sol_truth["wse_truth"]
    h_truth = sol_truth["h_truth"]
    u_truth = sol_truth["u_truth"]
    v_truth = sol_truth["v_truth"]
    zb_cells_truth = sol_truth["zb_cells_truth"]
    ManningN_values_truth = sol_truth["ManningN_values_truth"]
    inlet_discharges_truth = sol_truth["inlet_discharges_truth"]

    #combine the truth data into a dictionary
    observed_data = Dict("WSE_truth" => WSE_truth, "h_truth" => h_truth, "u_truth" => u_truth, "v_truth" => v_truth, "zb_cells_truth" => zb_cells_truth, "ManningN_values_truth" => ManningN_values_truth, "inlet_discharges_truth" => inlet_discharges_truth)


    #debug start
    # @show Q0
    # @show tspan
    # @show params_array
    # @assert !ismissing(params_array) "params_array contains missing values!"
    # @assert !isnothing(params_array) "params_array contains `nothing` values!"

    # prob = ODEProblem(ode_f, Q0, tspan, params_array)

    # #use Enzyme to test the gradient of the ODE and identify the source of the error
    # #See https://docs.sciml.ai/SciMLSensitivity/dev/faq/
    # SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true
    # p = prob.p
    # y = prob.u0
    # f = prob.f
    # t = tspan[1]  # Add this line to define t
    # @show typeof(p)
    # @show typeof(y)
    # @show typeof(f)
    # @show typeof(t)
    # @show p
    # @show y
    # @show f
    # @show t


    # # Test forward pass first
    # try
    #     #test_forward = f(y, p, t)
    #     #test_forward = swe_2d_rhs(y, p, active_range, param_ranges, t, settings,
    #     #    my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants,
    #     #    ManningN_cells, ManningN_ghostCells, inletQ_Length, inletQ_TotalQ, exitH_WSE,
    #     #    zb_cells, zb_ghostCells, zb_faces, S0)

    #     #test_forward = swe_2d_rhs(y, p, active_range, param_ranges, t, settings,
    #     #    my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants, inletQ_Length, exitH_WSE,
    #     #    )

    #     #@show typeof(test_forward)
    #     #@show size(test_forward)
    #     #@show test_forward

    #     println("Forward pass successful")
    #     println(" ")
    # catch e
    #     println("Forward pass failed")
    #     println(" ")
    #     @show e
    # end

    # #try
    #     # Test gradient computation directly
    #     #@code_warntype swe_2d_rhs(y, p, active_range, param_ranges, t, settings,
    #     #    my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants,
    #     #    ManningN_cells, ManningN_ghostCells, inletQ_Length, inletQ_TotalQ, exitH_WSE,
    #     #    zb_cells, zb_ghostCells, zb_faces, S0)
    #     SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true

    #     #use ForwardDiff to test the gradient of the ODE and identify the source of the error
    #     # gradient_output = ForwardDiff.gradient((y, p) -> begin
    #     #     #result = f(y, p, t)
    #     #     #result = swe_2d_ode(y, p, t)

    #     #     result = swe_2d_rhs(y, p, active_range, param_ranges, t, settings,
    #     #         my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants,
    #     #         ManningN_cells, ManningN_ghostCells, inletQ_Length, inletQ_TotalQ, exitH_WSE,
    #     #         zb_cells, zb_ghostCells, zb_faces, S0)

    #     #     @show typeof(result)
    #     #     @show size(result)
    #     #     @show result

    #     #     # Sum to get scalar output for gradient
    #     #     sum(result)
    #     # end, y, p)

    #     gradient_output = Zygote.gradient((y, p) -> begin
    #         #result = f(y, p, t)
    #         #result = swe_2d_ode(y, p, t)

    #         result = swe_2d_rhs(y, p, active_range, param_ranges, t, settings,
    #             my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants,
    #             ManningN_cells, ManningN_ghostCells, inletQ_Length, inletQ_TotalQ, exitH_WSE,
    #             zb_cells, zb_ghostCells, zb_faces, S0)

    #         #result = swe_2d_rhs(y, p, active_range, param_ranges, t, settings,
    #         #    my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants, inletQ_Length, exitH_WSE,
    #         #    )

    #         @show typeof(result)
    #         @show size(result)
    #         @show result

    #         # Sum to get scalar output for gradient
    #         sum(result)
    #         #return result[1,1]
    #     end, y, p)

    #     @show typeof(gradient_output)
    #     @show size.(gradient_output)
    #     @show gradient_output

    #     println("After Zygote.gradient, Gradient computation successful")
    #     println(" ")
    # #catch e
    # #    println("Gradient pass failed")

    # #    @show e

    # #    println(" ")
    # #end

    # # Now test the pullback with more detailed error catching
    # # try
    # #     λ = ones(size(prob.u0)) #zero(prob.u0)

    # #     #Zygote.pullback takes two arguments:
    # #     #  First argument: a function that we want to differentiate
    # #     #  Remaining arguments: the values at which to evaluate the function (y and p in this case)
    # #     #_dy is the result of the forward pass; back is the gradient function
    # #     #_dy, back = Zygote.pullback((u, p) -> f(u, p, t), y, p)
    # #     _dy, back = Zygote.pullback((u, p) -> Array(f(u, p, t)), y, p)

    # #     #_dy, back = Zygote.pullback(y, p) do u, p  
    # #     #    #vec(f(u, p, t))
    # #     #    f(u, p, t)
    # #     #end
    # #     println("Pullback creation successful")
    # #     @show typeof(_dy)
    # #     @show size(_dy)
    # #     @show _dy

    # #     try
    # #          # Convert λ to match _dy type
    # #         λ = convert(typeof(_dy), λ)

    # #         tmp1, tmp2 = back(λ)                  #tmp1 is the gradient of the state variables; tmp2 is the gradient of the parameters
    # #         println("Backward pass successful")
    # #         @show typeof(tmp1)
    # #         @show size(tmp1)
    # #         @show typeof(tmp2)
    # #         @show size(tmp2)
    # #         @show tmp1
    # #         @show tmp2
    # #     catch e
    # #         println("Backward pass failed")
    # #         @show e
    # #         @show typeof(λ)
    # #         @show size(λ)
    # #         @show λ
    # #     end
    # # catch e
    # #     println("Pullback creation failed")
    # #     @show e
    # # end

    # #throw("stop here")

    # #debug end
    # end



    # Define the loss function
    function compute_loss(p, Q0, tspan, observed_data, params_ranges, data_type)
        # Create ODEProblem
        #prob = ODEProblem(shallow_water_eq!, u0, tspan, p)
        prob = ODEProblem(ode_f, Q0, tspan, p)

        # Solve the ODE (forward pass)
        #For sensealg argument in ODE solve function: See https://docs.sciml.ai/SciMLSensitivity/dev/manual/differential_equation_sensitivities/ for the choice of AD type (sensealg)
        #If not specified, the default is a smart polyalgorithm used to automatically determine the most appropriate method for a given equation. It is defined in 
        #SciMLSensitivity.jl/src/concrete_solve.jl 
        #default_sensealg = if p !== SciMLBase.NullParameters() &&
        #                  !(eltype(u0) <: ForwardDiff.Dual) &&
        #                  !(eltype(p) <: ForwardDiff.Dual) &&
        #                  !(eltype(u0) <: Complex) &&
        #                  !(eltype(p) <: Complex) &&
        #                  length(u0) + length(tunables) <= 100
        #    ForwardDiffSensitivity()    #for small problems, it uses ForwardDiffSensitivity()
        #    ...
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
        if settings.inversion_settings.ode_solver == "Tsit5()"
            #pred = solve(prob, Tsit5(), adaptive=settings.inversion_settings.ode_solver_adaptive, dt=dt, saveat=t_save)  #working, but no control on sensealg
            #pred = solve(prob, Tsit5(), adaptive=settings.inversion_settings.ode_solver_adaptive, dt=dt, saveat=t_save, sensealg=ZygoteVJP())   #not working
            #pred = solve(prob, Tsit5(), adaptive=settings.inversion_settings.ode_solver_adaptive, dt=dt, saveat=t_save, sensealg=InterpolatingAdjoint(autojacvec=ZygoteVJP())) #working
            pred = solve(prob, Tsit5(), adaptive=settings.inversion_settings.ode_solver_adaptive, dt=dt, saveat=t_save, sensealg=ForwardDiffSensitivity()) #working
            #pred = solve(prob, Tsit5(), adaptive=settings.inversion_settings.ode_solver_adaptive, dt=dt, saveat=t_save, sensealg=ReverseDiffAdjoint()) #not working, ReverseDiffAdjoint only supports vector u0.
            #pred = solve(prob, Tsit5(), adaptive=settings.inversion_settings.ode_solver_adaptive, dt=dt, saveat=t_save, sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP())) #working
            #pred = solve(prob, Tsit5(), adaptive=settings.inversion_settings.ode_solver_adaptive, dt=dt, saveat=t_save, sensealg=ReverseDiffVJP()) #not working
        else
            #pred = solve(prob, Tsit5(), adaptive=true, dt=dt, saveat=t_save)
            error("Not implemented yet")
        end

        Zygote.ignore() do
            #@show pred.retcode
            #@show typeof(pred)
            #@show size(pred)
            #@show pred
        end

        #compute the loss
        # Ensure type stability in loss computation
        loss_total = zero(data_type)
        loss_pred = zero(data_type)
        loss_pred_WSE = zero(data_type)
        loss_pred_uv = zero(data_type)
        loss_slope = zero(data_type)

        if pred.retcode == SciMLBase.ReturnCode.Success
            WSE_truth = observed_data["WSE_truth"]
            h_truth = observed_data["h_truth"]
            u_truth = observed_data["u_truth"]
            v_truth = observed_data["v_truth"]

            #Get the bed elevation at cells
            zb_cells_temp = @view p[params_ranges.zb_start:params_ranges.zb_end]

            l = pred[:, 1, end] .+ zb_cells_temp .- WSE_truth  #loss = free surface elevation mismatch

            #loss for free surface elevation mismatch
            loss_pred_WSE = sum(abs2, l)

            #loss for velocity mismatch
            # Add small epsilon to prevent division by zero
            ϵ = sqrt(eps(data_type))

            if settings.inversion_settings.bInversion_u_loss      #if also include u in the loss 
                l_u = pred[:, 2, end] ./ (pred[:, 1, end] .+ ϵ) .- u_truth
                l_v = pred[:, 3, end] ./ (pred[:, 1, end] .+ ϵ) .- v_truth

                loss_pred_uv = sum(abs2, l_u) + sum(abs2, l_v)
            end

            #combined loss due to free surface elevation mismatch and velocity mismatch
            loss_pred = loss_pred_WSE + loss_pred_uv

            #loss for bed slope regularization
            if settings.inversion_settings.bInversion_slope_loss    #if bed slope is included in the loss 
                loss_slope = calc_slope_loss(zb_cells_temp, my_mesh_2D)
            end

            #combined loss due to free surface elevation mismatch, velocity mismatch, and bed slope regularization
            loss_total = loss_pred + loss_slope
        else
            loss_total = convert(data_type, Inf)
        end

        Zygote.ignore() do
            #@show loss_total
            #@show loss_pred
            #@show loss_pred_WSE
            #@show loss_pred_uv
            #@show loss_slope
        end

        return loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope, pred
    end

    function optimize_parameters(Q0, tspan, observed_data, p_init, active_range, param_ranges)
        # Loss function for optimization
        function opt_loss(θ, p)  # Add p argument even if unused

            data_type = eltype(θ)

            Zygote.ignore() do
                #println("data_type = ", data_type)
                #println("active_range = ", active_range)
                #println("param_ranges = ", param_ranges)
                #println("θ = ", θ)
                #println("p_init = ", p_init)
            end

            # Create new parameter set without mutation
            p_new = [i ∈ active_range ? θ[i - active_range[1] + 1] : p_init[i] for i in 1:length(p_init)]

            Zygote.ignore() do
                #println("p_new = ", p_new)
            end

            loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope, pred = compute_loss(p_new, Q0, tspan, observed_data, param_ranges, data_type)
            return loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope, pred
        end

        # Initial values for optimization parameters (vcat to flatten the array parameters to 1D array for optimizers)
        θ0 = p_init[active_range]  #get a copy of the subarray of p_init as the initial values for the optimization

        # Define AD type choice for optimization's gradient computation
        #The following is from SciMLSensitivity documentation regarding the choice of AD 
        #  AutoForwardDiff(): The fastest choice for small optimizations
        #  AutoReverseDiff(compile=false): A fast choice for large scalar optimizations
        #  AutoTracker(): Like ReverseDiff but GPU-compatible
        #  AutoZygote(): The fastest choice for non-mutating array-based (BLAS) functions
        #  AutoFiniteDiff(): Finite differencing, not optimal but always applicable
        #  AutoModelingToolkit(): The fastest choice for large scalar optimizations
        #  AutoEnzyme(): Highly performant AD choice for type stable and optimized code
        adtype = Optimization.AutoZygote()
        #adtype = Optimization.AutoReverseDiff(compile=false)
        #adtype = Optimization.AutoForwardDiff()

        # Define the optimization problem
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
        optf = OptimizationFunction(opt_loss, adtype)

        #From SciMLSensitivity documentation: https://docs.sciml.ai/Optimization/stable/API/optimization_problem/
        #OptimizationProblem{iip}(f, u0, p = SciMLBase.NullParameters(),;
        #                          lb = nothing,
        #                          ub = nothing,
        #                          lcons = nothing,
        #                          ucons = nothing,
        #                          sense = nothing,
        #                          kwargs...)
        optprob = OptimizationProblem(optf, θ0)  # No parameters needed

        # Define the bounds for the parameter (only applicable for some optimizers which support lb and ub)
        lb_p = zeros(my_mesh_2D.numOfCells)
        lb_p .= -0.1

        ub_p = zeros(my_mesh_2D.numOfCells)
        ub_p .= 0.3

        # Solve optimization problem
        # From SciMLSensitivity documentation: https://docs.sciml.ai/Optimization/stable/API/optimization_solution/
        # Returned optimization solution Fields:
        #   u: the representation of the optimization's solution.
        #   cache::AbstractOptimizationCache: the optimization cache` that was solved.
        #   alg: the algorithm type used by the solver.
        #   objective: Objective value of the solution
        #   retcode: the return code from the solver. Used to determine whether the solver solved successfully or whether 
        #            it exited due to an error. For more details, see the return code documentation.
        #   original: if the solver is wrapped from a external solver, e.g. Optim.jl, then this is the original return from said solver library.
        #   stats: statistics of the solver, such as the number of function evaluations required.
        if settings.inversion_settings.optimizer == "Adam"
            sol = solve(optprob, Adam(settings.inversion_settings.learning_rate), callback=callback, maxiters=settings.inversion_settings.max_iterations)
            #sol = solve(optprob, Adam(settings.inversion_settings.learning_rate), maxiters=settings.inversion_settings.max_iterations)
        elseif settings.inversion_settings.optimizer == "LBFGS"
            sol = solve(optprob, LBFGS(), callback=callback, maxiters=settings.inversion_settings.max_iterations)
        else
            error("Invalid optimizer choice. Supported optimizers: Adam, LBFGS. No inversion is performed.")
        end

        #timing and profiling tools
        #@time, #@profile and #Profile.print()
        #@show sol

        return sol
    end

    #define the accumulators for the inversion results
    LOSS = []                              # Loss accumulator
    PRED = []                              # prediction accumulator
    PARS = []                              # parameters accumulator

    callback = function (θ, loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope, pred) #callback function to observe training
        iter = size(LOSS)[1]  #get the inversion iteration number (=length of LOSS array)

        Zygote.ignore() do
            println("      iter, loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope = ", iter, ", ",
                loss_total, ", ", loss_pred, ", ", loss_pred_WSE, ", ", loss_pred_uv, ", ", loss_slope)
        end

        append!(PRED, [pred[:, 1, end]])
        append!(LOSS, [[loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope, pred]])

        if !isa(θ, Vector{Float64})  #NLopt returns an optimization object, not an arrary
            #println("theta.u = ", θ.u)
            append!(PARS, [copy(θ.u)])
        else
            append!(PARS, θ)
        end

        #if l > 1e-9
        #    false
        #else 
        #    true   #force the optimizer to stop 
        #end

        false
    end

    #inversion parameters: 
    # initial guess for the inversion parameters
    params_array_init = params_array

    #perform the inversion
    sol = optimize_parameters(Q0, tspan, observed_data, params_array_init, active_range, param_ranges)

    #save the inversion results
    jldsave(joinpath(save_path, settings.inversion_settings.save_file_name); LOSS, PRED, PARS)

    #process the inversion results
    println("   Post-processing inversion results ...")

    #process inversion results
    postprocess_inversion_results_2D(settings, my_mesh_2D, nodeCoordinates, zb_cells_truth, h_truth, u_truth, v_truth, WSE_truth)

end

#Timing 
end_time = now()  # Current date and time
elapsed_time = end_time - start_time
elapsed_seconds = Millisecond(elapsed_time).value / 1000
println("Elapsed time in seconds: $elapsed_seconds")

println("All done!")
