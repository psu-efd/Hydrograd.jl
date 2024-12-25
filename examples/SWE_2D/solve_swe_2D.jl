using Revise

using Dates

using JLD2

using IterTools

using PyCall

using Profile

using AdHydraulics

using OrdinaryDiffEq

using SparseArrays

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
using Enzyme

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
include("process_inversion_results_2D.jl")

println("Solving 2D SWE...")

#define control variables
bSimulate_Synthetic_Data = false    #whether to do the 1D SWE simulation to create synthetic data 
bPerform_Inversion = true           #whether to do inversion 
bPlot_Inversion_Results = true     #whehter to plot the inversion results

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

#create the 2D mesh 
my_mesh_2D = srh_all_Dict["my_mesh_2D"]

#mesh related variables which might be mutable (not in struct mesh_2D)
#nodeCoordinates::Array{Float64,3}  # Node coordinates: Float64 2D array [numOfNodes, 3]
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

#If simulating synthetic data, make a copy of the bathymetry truth: zb_cell_truth
if bSimulate_Synthetic_Data
    zb_cells_truth = deepcopy(zb_cells)
end

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
setup_initial_condition!(my_mesh_2D, nodeCoordinates, eta, zb_cells, h, q_x, q_y, true)

#setup ghost cells for initial condition
setup_ghost_cells_initial_condition!(my_mesh_2D, eta, h, q_x, q_y, eta_ghostCells, h_ghostCells, q_x_ghostCells, q_y_ghostCells)

#create and preprocess boundary conditions: boundary_conditions only contains the static information of the boundaries.
boundary_conditions, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length, inletQ_TotalA, inletQ_DryWet, 
exitH_WSE, exitH_H, exitH_A, wall_H, wall_A, symm_H, symm_A = initialize_boundary_conditions_2D(srh_all_Dict, nodeCoordinates)
#println("boundary_conditions = ", boundary_conditions)
#throw(ErrorException("stop here"))

#set up ODE parameters. In this problem, the "parameter" is para=zb_cells.
para = zb_cells
Q0 = hcat(h, q_x, q_y)   #initial condition: Q = [h q_x q_y]

Q_ghost = hcat(h_ghostCells, q_x_ghostCells, q_y_ghostCells)   #ghost cell values of Q

# time information 
#tspan = (swe_2D_constants.tStart, swe_2D_constants.tEnd)
#dt = swe_2D_constants.dt
#t = tspan[1]:dt:tspan[2]

tspan = (0.0, 100.0)    #100.0
dt = 0.1
t = tspan[1]:dt:tspan[2]

dt_save = (tspan[2] - tspan[1])/100.0
t_save = tspan[1]:dt_save:tspan[2]
println("t_save = ", t_save)

# Populate the Jacobian sparsity pattern
# Assume each variable depends on itself and its neighbors
# we have 3 variables (h, q_x, q_y) for each cell
jac_sparsity = spzeros(3*my_mesh_2D.numOfCells, 3*my_mesh_2D.numOfCells)

println("my_mesh_2D.cellNeighbors_Dict = ", my_mesh_2D.cellNeighbors_Dict)

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

#println("jac_sparsity = ", jac_sparsity)
#println("size of jac_sparsity = ", size(jac_sparsity))
#println("nnz of jac_sparsity = ", nnz(jac_sparsity))
#throw("stop here")

# Define the ODE function with explicit types
function swe_2d_ode(Q, para, t)
    dQdt = swe_2d_rhs(Q, para, t, my_mesh_2D, boundary_conditions, swe_2D_constants, ManningN_cells, inletQ_Length, inletQ_TotalQ, exitH_WSE)    

    return dQdt
end

# Create the ODEFunction with the typed function
#ode_f = ODEFunction(swe_2d_ode; jac_prototype=nothing)
ode_f = ODEFunction(swe_2d_ode; jac_prototype=jac_sparsity)

prob = ODEProblem(ode_f, Q0, tspan, para)

#use Enzyme to test the gradient of the ODE and identify the source of the error
#See https://docs.sciml.ai/SciMLSensitivity/dev/faq/
# SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true
# p = prob.p
# y = prob.u0
# f = prob.f
# λ = zero(prob.u0)
# _dy, back = Zygote.pullback(y, p) do u, p
#     vec(f(u, p, t))
# end
# tmp1, tmp2 = back(λ)

# @show tmp1
# @show tmp2

# throw("stop here")

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
        jldsave(joinpath(save_path, "zb_cells_truth.jld2"); zb_cells_truth)

        swe_2D_save_results_SciML(sol, total_water_volume, my_mesh_2D, nodeCoordinates, zb_cells, save_path)
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

    zb_cells_truth = load(joinpath(save_path, "zb_cells_truth.jld2"))["zb_cells_truth"]

    WSE_truth = h_truth .+ zb_cells_truth

    #inversion parameters: zb_cell
    # initial guess for the parameters
    #ps = zeros(Float64, my_mesh_2D.numOfCells) .+ 0.01
    ps = zeros(Float64, my_mesh_2D.numOfCells) 
    
    function predict(θ)
        #Zygote.ignore() do
        #    println("Current parameters: ", ForwardDiff.value.(θ))
        #end

        #See https://docs.sciml.ai/SciMLSensitivity/dev/faq/ for the choice of AD type (sensealg)
        #If not specified, the default is a smart polyalgorithm is used to automatically determine the most appropriate method for a given equation.
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
        #sol = solve(prob, Tsit5(), adaptive=false, p=θ, dt=dt, saveat=t_save)  #[:,1,end]  #not working for long time span
        #sol = solve(prob, Euler(), adaptive=false, p=θ, dt=dt, saveat=t_save)  #[:,1,end]

        #sol = solve(prob, Tsit5(), p=θ, dt=dt, saveat=t_save; sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP())) #only works for short time span
        #sol = solve(prob, Tsit5(), p=θ, dt=dt, saveat=t_save; sensealg=ForwardDiffSensitivity())   #runs, but very slow
        #sol = solve(prob, Tsit5(), adaptive=false, p=θ, dt=dt, saveat=t_save; sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()))
        sol = solve(prob, Tsit5(), adaptive=false, p=θ, dt=dt, saveat=t_save; sensealg=GaussAdjoint(autojacvec=ZygoteVJP()))

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
        #l = sol[:,1,end] .- 0.5  #loss = free surface elevation mismatch
        l = sol[:,1,end] .- h_truth
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

    #function test_ode_gradient(θ)
    #    sol = solve(prob, Tsit5(), p=θ, dt=dt, saveat=t_save; sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()))
    #    h_pred = Array(sol)[:,1,end]
    #    return sum(abs2, h_pred)
    #end
    
    #grads = Zygote.gradient(test_ode_gradient, ps)
    #println(grads)

    #throw("stop here")

    # Add this before gradient computation
    #@code_warntype(my_loss(ps))

    SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true
    #grad = ForwardDiff.gradient(my_loss, ps)
    #grad = Zygote.gradient(my_loss, ps)
    #tape = ReverseDiff.GradientTape(my_loss, ps)
    #ReverseDiff.compile!(tape)  # Optionally compile the tape for better performance
    #grad = zeros(length(ps))  # Preallocate gradient array
    #grad = ReverseDiff.gradient(my_loss, ps)

    #jac = Zygote.jacobian(predict, ps)
    #jac = ReverseDiff.jacobian(predict, ps)
    #@show grad
    #plot(jac)
    #readline()
    #exit()
    #throw("stop here")

    ## Defining Loss function
    #zb_cell_local = zeros(Number, mesh.nCells)

    function loss(θ)
        pred = predict(θ)            #Forward prediction with current θ (=zb at cells)
        #l = pred[:,1,end] - 0.5 #h_truth  #loss = free surface elevation mismatch
        
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

    #grad = Zygote.gradient(loss, ps)
    #grad = ForwardDiff.gradient(loss, ps)
    #grad = ForwardDiff.gradient(predict, θ)
    #@show grad
    #println("grad = ", grad)
    #readline()
    #throw("stop here")


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

    #setup the bounds for the bathymetry
    lb_p = zeros(my_mesh_2D.numOfCells)
    lb_p .= -0.1  

    ub_p = zeros(my_mesh_2D.numOfCells)
    ub_p .= 0.3 

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
    #@show PARS 
    jldsave(joinpath(save_path, "inversion_results.jld2"); LOSS, PRED, PARS)

end

if bPlot_Inversion_Results

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
