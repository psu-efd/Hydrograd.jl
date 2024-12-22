using Revise

using Dates

using JLD2

using IterTools

using PyCall

using AdHydraulics

using OrdinaryDiffEq

#SciML
using SciMLSensitivity

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

#for Bayesian estimation
#using Turing

# Load StatsPlots for visualizations and diagnostics.
#using StatsPlots

#using LinearAlgebra

using Plots

using InteractiveUtils

#using Random
#Random.seed!(1234)

#Timing 
start_time = now()  # Current date and time

include("process_SRH_2D_input.jl")
include("process_ManningN_2D.jl")
include("process_BCs_2D.jl")
include("process_bed_2D.jl")
include("process_ICs_2D.jl")
include("semi_discretize_2D.jl")
include("misc_utilities_2D.jl")
include("custom_ODE_solver.jl")
println("Solving 2D SWE...")

#define control variables
bSimulate_Synthetic_Data = true    #whether to do the 1D SWE simulation to create synthetic data 
bPlot_Simulation_Results = false   #whether to plot simulation results
bPerform_Inversion = false           #whether to do inversion 
bPlot_Inversion_Results = false     #whehter to plot the inversion results

#options for inversion 
bInversion_slope_loss = false   #whether to include slope loss 
bInversion_include_u = false    #whehter to add velocity to the loss function (by default, we already have water surface elevatino eta)

#directory to save to (the same directory as the main file)
save_path = dirname(@__FILE__)

#define a swe_2D_constants object with some default values
swe_2D_constants = swe_2D_consts(t=0.0, dt=0.1, tStart=0.0, tEnd=1.0)

#read data from SRH-2D hydro, geom, and material files
#srhhydro_file_name = "simple.srhhydro"
srhhydro_file_name = "oneD_channel_with_bump.srhhydro"
#srhhydro_file_name = "oneD_channel_with_bump_all_walls.srhhydro"
#srhhydro_file_name = "twoD_channel_with_bump.srhhydro"

srh_all_Dict = process_SRH_2D_input(srhhydro_file_name)

#update swe_2D_constants based on the SRH-2D data
update_swe_2D_constants!(swe_2D_constants, srh_all_Dict)

#get the 2D mesh 
my_mesh_2D = srh_all_Dict["my_mesh_2D"]

#  setup bed elevation 
#zb_faces = zeros(Float64, my_mesh_2D.numOfFaces)      #zb at faces 
#zb_cells = zeros(Float64, my_mesh_2D.numOfCells)      #zb at cell centers 
#zb_ghostCells = zeros(Float64, my_mesh_2D.numOfAllBounaryFaces)   #zb at ghost cell centers 
#S0 = zeros(Float64, my_mesh_2D.numOfCells, 2)          #bed slope at cell centers 

#If performing inversion, set the bed elevation to zero
if bPerform_Inversion
    my_mesh_2D.nodeCoordinates[:,3] .= 0.0
end

#setup bed elevation: computer zb at cell centers from nodes, then interpolate zb from cell to face and compute the bed slope at cells
zb_cells, zb_ghostCells, zb_faces, S0 = setup_bed(my_mesh_2D, true)

#make a copy of the bathymetry truth: zb_cell_truth
zb_cells_truth = deepcopy(zb_cells)

#setup Manning's n 
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
setup_initial_condition!(my_mesh_2D, eta, zb_cells, h, q_x, q_y, true)

#setup ghost cells for initial condition
setup_ghost_cells_initial_condition!(my_mesh_2D, eta, h, q_x, q_y, eta_ghostCells, h_ghostCells, q_x_ghostCells, q_y_ghostCells)

#preprocess boundary conditions 
inletQ_BC_indices = Int[]   #indices of inlet-Q boundaries in the boundary list
exitH_BC_indices = Int[]    #indices of exit-H boundaries in the boundary list
wall_BC_indices = Int[]     #indices of wall boundaries in the boundary list
symm_BC_indices = Int[]     #indices of symmetry (no slip) boundaries in the boundary list

#compute the indices of each boundary in the global boundary list (nWall_BCs and nSymm_BCs are updated in this function)
compute_boundary_indices!(my_mesh_2D, srh_all_Dict, inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices)

nInletQ_BCs = srh_all_Dict["nInletQ_BCs"]   #number of inlet-Q boundaries
nExitH_BCs = srh_all_Dict["nExitH_BCs"]     #number of exit-H boundaries
nWall_BCs = srh_all_Dict["nWall_BCs"]       #number of wall boundaries
nSymm_BCs = srh_all_Dict["nSymm_BCs"]       #number of symmetry boundaries

#some data arrays for inlet-q bounaries:
inletQ_faceIDs = Array{Array{Int}}(undef, nInletQ_BCs)   #face IDs of the inlet-q boundaries
inletQ_ghostCellIDs = Array{Array{Int}}(undef, nInletQ_BCs)   #ghost cell IDs of the inlet-q boundaries
inletQ_internalCellIDs = Array{Array{Int}}(undef, nInletQ_BCs)   #internal cell IDs of the inlet-q boundaries

inletQ_faceCentroids = Vector{Matrix{Float64}}(undef, nInletQ_BCs)   #face centroids of the inlet-q boundaries
inletQ_faceOutwardNormals = Vector{Matrix{Float64}}(undef, nInletQ_BCs)   #face outward normals of the inlet-q boundaries

inletQ_TotalQ = zeros(Float64, nInletQ_BCs)   #total discharge for each inlet-q boundary
inletQ_H = Array{Array{Float64}}(undef, nInletQ_BCs)   #inlet water depth for each face in each inlet-q boundary
inletQ_A = Array{Array{Float64}}(undef, nInletQ_BCs)   #area for each face in each inlet-q boundary
inletQ_ManningN = Array{Array{Float64}}(undef, nInletQ_BCs)   #Manning's n for each face in each inlet-q boundary
inletQ_Length = Array{Array{Float64}}(undef, nInletQ_BCs)   #length for each face in each inlet-q boundary
inletQ_TotalA = zeros(Float64, nInletQ_BCs)   #total cross-sectional area for each inlet-q boundary
inletQ_DryWet = Array{Array{Int}}(undef, nInletQ_BCs)   #dry(=0)/wet(=1) flag for each inlet-q boundary

#preprocess inlet-q boundaries
preprocess_inlet_q_boundaries(my_mesh_2D, nInletQ_BCs, srh_all_Dict["srhhydro_IQParams"], inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
    inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, 
    inletQ_ManningN, inletQ_Length, inletQ_TotalA, inletQ_DryWet)

#some data arrays for exit-h bounaries:
exitH_faceIDs = Array{Array{Int}}(undef, nExitH_BCs)   #face IDs of the exit-h boundaries
exitH_ghostCellIDs = Array{Array{Int}}(undef, nExitH_BCs)   #ghost cell IDs of the exit-h boundaries
exitH_internalCellIDs = Array{Array{Int}}(undef, nExitH_BCs)   #internal cell IDs of the exit-h boundaries

exitH_faceCentroids = Vector{Matrix{Float64}}(undef, nExitH_BCs)   #face centroids of the exit-h boundaries
exitH_faceOutwardNormals = Vector{Matrix{Float64}}(undef, nExitH_BCs)   #face outward normals of the exit-h boundaries

exitH_WSE = Array{Float64}(undef, nExitH_BCs)   #WSE for each exit-h boundary

exitH_H = Array{Array{Float64}}(undef, nExitH_BCs)   #inlet water depth for each exit-h boundary
exitH_A = Array{Array{Float64}}(undef, nExitH_BCs)   #inlet cross-sectional area for each exit-h boundary

#preprocess exit-h boundaries
preprocess_exit_h_boundaries(my_mesh_2D, nExitH_BCs, srh_all_Dict["srhhydro_EWSParamsC"], exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
    exitH_internalCellIDs, exitH_faceCentroids, exitH_faceOutwardNormals, exitH_WSE, exitH_H, exitH_A)

#some data arrays for wall bounaries:
wall_faceIDs = Array{Array{Int}}(undef, nWall_BCs)   #face IDs of the wall boundaries
wall_ghostCellIDs = Array{Array{Int}}(undef, nWall_BCs)   #ghost cell IDs of the wall boundaries
wall_internalCellIDs = Array{Array{Int}}(undef, nWall_BCs)   #internal cell IDs of the wall boundaries

wall_faceCentroids = Vector{Matrix{Float64}}(undef, nWall_BCs)   #face centroids of the wall boundaries
wall_outwardNormals = Vector{Matrix{Float64}}(undef, nWall_BCs)   #face outward normals of the wall boundaries

wall_H = Array{Array{Float64}}(undef, nWall_BCs)   #water depth for each wall boundary
wall_A = Array{Array{Float64}}(undef, nWall_BCs)   #cross-sectional area for each wall boundary

#index array for each wall boundary: if the index is in the wall boundary, the value is 1.0, otherwise 0.0
#wall_idx = Array{Array{Float64}}(undef, nWall_BCs)   

#preprocess wall boundaries
preprocess_wall_boundaries(my_mesh_2D, nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs,
    wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals, wall_H, wall_A)

#some data arrays for symmetry boundaries:
symm_faceIDs = Array{Array{Int}}(undef, nSymm_BCs)
symm_ghostCellIDs = Array{Array{Int}}(undef, nSymm_BCs)
symm_internalCellIDs = Array{Array{Int}}(undef, nSymm_BCs)

symm_faceCentroids = Vector{Matrix{Float64}}(undef, nSymm_BCs)
symm_outwardNormals = Vector{Matrix{Float64}}(undef, nSymm_BCs)

symm_H = Array{Array{Float64}}(undef, nSymm_BCs)
symm_A = Array{Array{Float64}}(undef, nSymm_BCs)

#preprocess symmetry boundaries
preprocess_symmetry_boundaries(my_mesh_2D, nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs,
    symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals, symm_H, symm_A)

#println("inletQ_faceIDs: ", inletQ_faceIDs)
#println("exitH_faceIDs: ", exitH_faceIDs)
#println("wall_faceIDs: ", wall_faceIDs)
#println("symm_faceIDs: ", symm_faceIDs)

#set up ODE parameters. In this problem, the "parameter" is para=zb_cells.
para = zb_cells
Q0 = hcat(h, q_x, q_y)   #initial condition: Q = [h q_x q_y]

Q_ghost = hcat(h_ghostCells, q_x_ghostCells, q_y_ghostCells)   #ghost cell values of Q

# time information 
#tspan = (swe_2D_constants.tStart, swe_2D_constants.tEnd)
#dt = swe_2D_constants.dt
#t = tspan[1]:dt:tspan[2]

tspan = (0.0, 1.0)    #100.0
dt = 0.1
t = tspan[1]:dt:tspan[2]

dt_save = (tspan[2] - tspan[1])/10.0
t_save = tspan[1]:dt_save:tspan[2]
println("t_save = ", t_save)

# Define the ODE function with explicit types
function swe_2d_ode(Q, para, t)
    dQdt = swe_2d_rhs(Q, Q_ghost, para, t, my_mesh_2D, zb_cells, zb_ghostCells, zb_faces, S0,
            swe_2D_constants, ManningN_cells, ManningN_ghostCells,
            nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
            inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, 
            inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length,
            inletQ_TotalA, inletQ_DryWet,  
            nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
            exitH_internalCellIDs, exitH_faceCentroids, exitH_WSE,
            nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, 
            wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals,
            nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs, 
            symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals)    

    return dQdt
end

# Create the ODEFunction with the typed function
ode_f = ODEFunction(swe_2d_ode;
    jac_prototype=nothing)

prob = ODEProblem(ode_f, Q0, tspan, para)

if bSimulate_Synthetic_Data

    solver_choice = "SciML"
    #solver_choice = "MyOwn"

    println("   Performing 2D SWE simulation ...")

    #for direct sensivity analysis
    #prob = ODEForwardSensitivityProblem((Q, p, t) -> swe_1D_rhs!(Q, p, t, 
    #                          mesh, swe_1d_constants, left_bcType, right_bcType, S0), Q0, tspan, p)

    # solve the ODE with a solver of your choice 
    #bm = @benchmark solve(prob, Heun(), adaptive=false, dt=dt, saveat=t)
    #@show bm 

    if solver_choice == "SciML"

        println("   Performing 2D SWE simulation with SciML solver ...")
    
        sol = solve(prob, Tsit5(), adaptive=false, dt=dt, saveat=t_save)
        #sol = solve(prob, Heun(), adaptive=false, dt=dt, saveat=t)
        #sol = solve(prob, AB3(), adaptive=false, dt=dt, saveat = t)
        #sol = solve(prob, Rosenbrock23(), adaptive=false, dt=dt, saveat = t)  #implicit

        #x, dp = extract_local_sensitivities(sol)
        #da = dp[1]
        #plot(sol.t, da', lw = 3)

        #println("Press any key to exit ...")
        #readline()
        #exit()

        # #save the results
        # #save the simulation solution results
        jldsave(joinpath(save_path, "simulation_solution.jld2"); sol)

        swe_2D_save_results_SciML(sol, total_water_volume, my_mesh_2D, zb_cells, save_path)
    end

    #My own ODE Solver
    if solver_choice == "MyOwn"

        println("   Performing 2D SWE simulation with MyOwn solver ...")

        
        sol = my_solve(para, Q0, my_mesh_2D, tspan, dt)

        # #save the results
        # #save the simulation solution results
        jldsave(joinpath(save_path, "simulation_solution.jld2"); sol)

        swe_2D_save_results_custom(sol, total_water_volume, my_mesh_2D, zb_cells, save_path)
                
    end

    

end

if bPlot_Simulation_Results

    println("   Plotting 2D SWE simulation results ...")

    #open the simulation result
    sol = load(joinpath(save_path, "simulation_solution.jld2"))["sol"]

    #plot the results at tEnd
    #swe_1D_make_plots(save_path)

    #make an animation of the simulation as a function of time 
    #swe_1D_make_animation(sol, mesh, zb_cell, save_path)
end


################
#Inversion part# 
################

if bPerform_Inversion

    println("   Performing inversion ...")

    #open the simulation result (as the ground truth)
    sol = load(joinpath(save_path, "simulation_solution.jld2"))["sol"]
    h_truth = Array(sol)[:,1,end]
    u_truth = Array(sol)[:,2,end]./Array(sol)[:,1,end]
    v_truth = Array(sol)[:,3,end]./Array(sol)[:,1,end]
    #@show sol

    #inversion parameters: zb_cell
    # initial guess for the parameters
    ps = zeros(Float64, my_mesh_2D.numOfCells) .+ 0.01
    
    function predict(θ)
        #Array(solve(prob, Heun(), adaptive=false, p=θ, dt=dt, saveat=t))[:,1,end]
        sol = solve(prob, Tsit5(), adaptive=false, p=θ, dt=dt, saveat=t_save)  #[:,1,end]
        #sol = solve(prob, Euler(), adaptive=false, p=θ, dt=dt, saveat=t_save)  #[:,1,end]

        #solve the ODE with my own solver
        #sol = my_solve(θ, Q0, my_mesh_2D, tspan, dt)

        #if !SciMLBase.successful_retcode(sol)
        #    # Return a high cost instead of NaN
        #    return fill(convert(eltype(θ), Inf), size(sol[1]))
        #end
        
        return sol
    end

    #my_loss for debugging
    function my_loss(θ)
        # Add debug prints
        #Zygote.ignore() do
            #println("Current parameters: ", ForwardDiff.value.(θ))
        #end

        sol = Array(predict(θ))            #Forward prediction with current θ (=zb at cells)
        l = sol[:,1,end] .- 0.5  #loss = free surface elevation mismatch
        loss = sum(abs2, l)

        #Zygote.ignore() do
        #    println("loss value = ", ForwardDiff.value(loss))
        #    println("dloss/dθ = ", ForwardDiff.partials(loss))
        #    #println("loss = ", loss)
        #end

        #pred_real = ForwardDiff.value.(sol)
        #pred_real = sol
        
        #println("size of pred_real", size(pred_real))
        #println("size of h_truth", size(h_truth))
        #println("pred = ", pred_real[:,1,end])
        #println("h_truth = ", h_truth)

        #Zygote.ignore() do
            # plt = plot(pred_real[:,1,end], 
            #     label="Predicted",
            #     #ylim=(0, 10),
            #     xlabel="x",
            #     ylabel="h",
            #     linestyle=:solid)

            # plot!(h_truth, label="True", linestyle=:dash)

            # display(plt)
        #end

        return loss
    end

    #test predict
    #sol_value = predict(ps)
    #@show sol_value
    #println("sol_value = ", sol_value)

    #loss_value = my_loss(ps)
    #println("loss_value = ", loss_value)

    # Add this before gradient computation
    #@code_warntype(my_loss(ps))

    #SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true
    #grad = ForwardDiff.gradient(my_loss, ps)
    grad = Zygote.gradient(my_loss, ps)[1]
    #tape = ReverseDiff.GradientTape(my_loss, ps)
    #ReverseDiff.compile!(tape)  # Optionally compile the tape for better performance
    #grad = zeros(length(ps))  # Preallocate gradient array
    #grad = ReverseDiff.gradient(my_loss, ps)

    #jac = Zygote.jacobian(predict, ps)
    #jac = ReverseDiff.jacobian(predict, ps)
    @show grad
    #plot(jac)
    #readline()
    #exit()
    throw("stop here")

    ## Defining Loss function
    #zb_cell_local = zeros(Number, mesh.nCells)

    function loss(θ)
        pred = predict(θ)            #Forward prediction with current θ (=zb at cells)
        l = pred[:,1,end] - 0.5 #h_truth  #loss = free surface elevation mismatch
        loss_pred_eta = sum(abs2, l)

        # loss_pred_uv = zero(eltype(θ))

        #  # Add small epsilon to prevent division by zero
        # ϵ = sqrt(eps(eltype(θ)))
        
        # if bInversion_include_u      #if also include u in the loss 
        #     l_u = pred[:,2,end]./(pred[:,1,end] .+ ϵ) .- u_truth
        #     l_v = pred[:,3,end]./(pred[:,1,end] .+ ϵ) .- v_truth

        #     loss_pred_uv = sum(abs2, l_u) + sum(abs2, l_v)
        # end 

        # loss_pred = loss_pred_eta + loss_pred_uv

        # loss_slope = zero(eltype(θ))

        # if bInversion_slope_loss    #if bed slope is included in the loss 
        #     loss_slope = calc_slope_loss(θ, my_mesh_2D)
        # end 

        # loss_total = loss_pred + loss_slope

        # #return loss_total, loss_pred, loss_pred_eta, loss_pred_uv, loss_slope, pred
        # return loss_total

        Zygote.ignore() do
            println("loss_pred_eta = ", loss_pred_eta)
        end

        return loss_pred_eta
    end

    grad = Zygote.gradient(loss, ps)
    #grad = ForwardDiff.gradient(loss, ps)
    #grad = ForwardDiff.gradient(predict, θ)
    @show grad
    #println("grad = ", grad)
    #readline()
    throw("stop here")


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

    # Let see prediction vs. Truth
    #scatter(sol[:, end], label = "Truth", size = (800, 500))
    #plot!(PRED[end][:, end], lw = 2, label = "Prediction")

    #The following is from SciMLSensitivity documentation regarding the choice of AD 
    #  AutoForwardDiff(): The fastest choice for small optimizations
    #  AutoReverseDiff(compile=false): A fast choice for large scalar optimizations
    #  AutoTracker(): Like ReverseDiff but GPU-compatible
    #  AutoZygote(): The fastest choice for non-mutating array-based (BLAS) functions
    #  AutoFiniteDiff(): Finite differencing, not optimal but always applicable
    #  AutoModelingToolkit(): The fastest choice for large scalar optimizations
    #  AutoEnzyme(): Highly performant AD choice for type stable and optimized code

    adtype = Optimization.AutoZygote()         #works on MBP, but slow. Not on Windows (don't know why). 
    
    #adtype = Optimization.AutoForwardDiff()  #ForwardDiff.jl is not supported for this problem.
    #adtype = Optimization.AutoReverseDiff()
    #adtype = Optimization.AutoEnzyme()   #failed, don't use "not implemented yet error".

    optf = Optimization.OptimizationFunction((θ, p) -> loss(θ), adtype)

    #setup the bounds for the bathymetry
    lb_p = zeros(my_mesh_2D.numOfCells)
    lb_p .= -0.1  

    ub_p = zeros(my_mesh_2D.numOfCells)
    ub_p .= 0.3 

    #optprob = Optimization.OptimizationProblem(optf, ps, lb=lb_p, ub=ub_p)
    optprob = Optimization.OptimizationProblem(optf, ps)

    #res = Optimization.solve(optprob, PolyOpt(), callback = callback)  #PolyOpt does not support lb and ub 
    #res = Optimization.solve(optprob, NLopt.LD_LBFGS(), callback = callback)   #very fast 
    #res = Optimization.solve(optprob, Optim.BFGS(), callback=callback; iterations=30, maxiters=40, f_calls_limit=20, show_trace=true)
    #res = Optimization.solve(optprob, Optim.BFGS(), callback=callback, maxiters = 100; show_trace=false)  #f_tol=1e-3, iterations=10, local_maxiters = 10
    #res = Optimization.solve(optprob, Optim.LBFGS(), callback=callback)  #oscilates around 1e-7
    #res = Optimization.solve(optprob, Optim.Newton(), callback=callback)  #error: not supported as the Fminbox optimizer
    #res = Optimization.solve(optprob, Optim.GradientDescent(), callback=callback)  #very slow decrease in loss 
    res = Optimization.solve(optprob, Adam(0.01), callback=callback, maxiters=100)
    
    @show res
    @show res.original

    #save the inversion results
    #@show PARS 
    jldsave(joinpath(save_path, "inversion_results.jld2"); LOSS, PRED, PARS, sol)

end

#Timing 
end_time = now()  # Current date and time
elapsed_time = end_time - start_time
elapsed_seconds = Millisecond(elapsed_time).value / 1000
println("Elapsed time in seconds: $elapsed_seconds")

println("All done!")
