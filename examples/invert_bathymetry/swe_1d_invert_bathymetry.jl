#This code solvers a 1D SWE with the following:
# BCs: inletQ on the left and exitH on the right
# Bed profile: prescribed 
# ICs: prescribed 
# 
#Then the code does an inversion of the bathymetry

using Dates
using Plots
using JLD2
using BenchmarkTools

using AdHydraulics

#SciML related packages
# DE
using DifferentialEquations
using OrdinaryDiffEq

#SciML
using SciMLSensitivity

#Optimizers
using Optimization
using OptimizationPolyalgorithms
using OptimizationNLopt
using Optim
using OptimizationFlux

#AD engines
using Zygote
using ForwardDiff
using ReverseDiff

#for Bayesian estimation
#using Turing

# Load StatsPlots for visualizations and diagnostics.
#using StatsPlots

#using LinearAlgebra

#using Random
#Random.seed!(1234)

#include problem-specific Julia files (must be in the same directory as the main file)
include("process_bed.jl")
include("process_initial_condition.jl")
include("semi_discretize.jl")
include("plot_inversion_results.jl")


time_start = Dates.now()

#print some information
println("SWE-1D simulation and inversion ...")

#directory to save to (the same directory as the main file)
save_path = dirname(@__FILE__)

#define control variables
bSimulate_Synthetic_Data = true    #whether to do the 1D SWE simulation to create synthetic data 
bPlot_Simulation_Results = true    #whether to plot simulation results
bPerform_Inversion = true           #whether to do inversion 
bPlot_Inversion_Results = true     #whehter to plot the inversion results

#options for inversion 
bInversion_slope_loss = true   #whether to include slope loss 
bInversion_include_u = true    #whehter to add velocity to the loss function (by default, we already have water surface elevatino eta)

#mesh related variables
nCells = 30     # Number of cells
L = 10.0  # Length of the domain

#create the 1D mesh
mesh = initialize_mesh_1D(nCells, L)

#create problem constants 
swe_1d_constants = swe_1D_const(
    g=9.81,
    t=0.0,             #starting time (current time)
    dt_min=0.001,      #minimum dt (if adaptive time step)
    dt=0.005,          #time step size
    CFL=0.4,           #CFL number 
    tEnd=0.1,          #end of simulation time 
    h_small=1.0e-3,      #a small water depth for dry/wet treatment
    RiemannSolver="HLL"  #choice of Riemann problem solver 
)

#variables

ManningN = 0.03    #Manning's n

#  setup bed 
#zb_face = zeros(mesh.nFaces)      #zb at faces (points in 1D)
zb_cell = zeros(Number, mesh.nCells)      #zb at cell centers 
#S0 = zeros(mesh.nCells)          #bed slope at cell centers 

setup_bed!(mesh, zb_cell)

#make a copy of the bathymetry truth: zb_cell_truth
zb_cell_truth = deepcopy(zb_cell)

# setup initial conditions 
eta = zeros(mesh.nCells)          #free surface elevation at cells 
h = zeros(mesh.nCells)            #water depth at cells 
q = zeros(mesh.nCells)            #q=hu at cells 
total_water_volume = Float64[]      #total volume of water in the domain 

setup_initial_eta!(mesh, eta, zb_cell, h, swe_1d_constants)

#setup boundary conditions
left_bcType = "inletQ"
#left_bcType = "wall"
left_bcValue = inletQ = 0.8
right_bcType = "exitH"
#right_bcType = "wall"
right_bcValue = exitH = 0.75

#fluxes on faces 
#fluxes = zeros(Number, 2, mesh.nFaces)
#temp flux
#flux = zeros(Number, 2)

#set up ODE parameters. In this problem, p=zb_cell
p = Float64.(zb_cell)
Q0 = hcat(h, q)   #initial condition: Q = [h q]

# time information 
tspan = (0.0, swe_1d_constants.tEnd)
dt = swe_1d_constants.dt
t = 0.0:dt:swe_1d_constants.tEnd

dt_save = (swe_1d_constants.tEnd - 0.0)/10.0
t_save = 0.0:dt_save:swe_1d_constants.tEnd

# define the ODE
f = ODEFunction((dQdt, Q, p, t) -> swe_1D_rhs!(dQdt, Q, p, t,
        mesh, swe_1d_constants, left_bcType, right_bcType, left_bcValue, right_bcValue, ManningN), 
        jac_prototype=nothing)

prob = ODEProblem(f, Q0, tspan, p)

if bSimulate_Synthetic_Data

    println("   Performing 1D SWE simulation ...")

    #for direct sensivity analysis
    #prob = ODEForwardSensitivityProblem((Q, p, t) -> swe_1D_rhs!(Q, p, t, 
    #                          mesh, swe_1d_constants, left_bcType, right_bcType, S0), Q0, tspan, p)

    # solve the ODE with a solver of your choice 
    #bm = @benchmark solve(prob, Heun(), adaptive=false, dt=dt, saveat=t)
    #@show bm 
    
    sol = solve(prob, Tsit5(), adaptive=true, dt=dt, saveat=t)
    #sol = solve(prob, Heun(), adaptive=false, dt=dt, saveat=t)
    #sol = solve(prob, AB3(), adaptive=false, dt=dt, saveat = t)
    #sol = solve(prob, Rosenbrock23(), adaptive=false, dt=dt, saveat = t)  #implicit

    #x, dp = extract_local_sensitivities(sol)
    #da = dp[1]
    #plot(sol.t, da', lw = 3)

    #println("Press any key to exit ...")
    #readline()
    #exit()

    #save the simulation solution results
    jldsave(joinpath(save_path, "simulation_solution.jld2"); sol)

    swe_1D_save_results(sol, total_water_volume, mesh, zb_cell, save_path)
end

if bPlot_Simulation_Results

    println("   Plotting 1D SWE simulation results ...")

    #open the simulation result
    sol = load(joinpath(save_path, "simulation_solution.jld2"))["sol"]

    #plot the results at tEnd
    swe_1D_make_plots(save_path)

    #make an animation of the simulation as a function of time 
    #swe_1D_make_animation(sol, mesh, zb_cell, save_path)
end

#println("Press enter to exit ...")
#readline()
#exit()

################
#Inversion part# 
################

if bPerform_Inversion

    println("   Performing inversion ...")

    #open the simulation result (as the ground truth)
    sol = load(joinpath(save_path, "simulation_solution.jld2"))["sol"]
    #@show sol

    #inversion parameters: zb_cell
    # initial guess for the parameters
    ps = zeros(mesh.nCells)

    function predict(θ)
        #Array(solve(prob, Heun(), adaptive=false, p=θ, dt=dt, saveat=t))[:,1,end]
        Array(solve(prob, Tsit5(), adaptive=true, p=θ, dt=dt, saveat=t))  #[:,1,end]
    end

    SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true
    #jac = ForwardDiff.jacobian(predict, ps)
    #jac = Zygote.jacobian(predict, ps)
    #jac = ReverseDiff.jacobian(predict, ps)
    #@show jac
    #plot(jac)
    #exit()

    ## Defining Loss function
    #zb_cell_local = zeros(Number, mesh.nCells)

    function loss(θ)
        pred = predict(θ)            #Forward prediction with current θ (=zb at cells)
        l = pred[:,1,end] - Array(sol)[:,1,end] + (θ - zb_cell_truth)  #loss = free surface elevation mismatch

        loss_pred = sum(abs2, l)
        loss_slope = calc_slope_loss(θ, mesh.dx)
        #loss_total = loss_pred + loss_slope
        loss_total = loss_pred

        #return loss_total, loss_pred, loss_slope, pred # Mean squared error + slope loss 
        return loss_pred, loss_pred, loss_slope, pred
    end

    #grad = Zygote.gradient(loss, ps)
    #grad = ForwardDiff.gradient(loss, ps)
    #grad = ForwardDiff.gradient(predict, θ)
    #@show grad
    #println("grad = ", grad)
    #readline()
    #exit()


    LOSS = []                              # Loss accumulator
    PRED = []                              # prediction accumulator
    PARS = []                              # parameters accumulator

    callback = function (θ, loss_total, loss_pred, loss_slope, pred) #callback function to observe training
        iter = size(LOSS)[1]  #get the inversion iteration number (=length of LOSS array)
        println("      iter, loss_total, loss_pred, loss_slope = ", iter, ", ", loss_total, ", ", loss_pred, ", ", loss_slope)

        append!(PRED, [pred[:,1,end]])
        append!(LOSS, [[loss_total, loss_pred, loss_slope]])

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

    adtype = Optimization.AutoZygote()         #works, but slow. Maybe because of memory allocations?
    #adtype = Optimization.AutoForwardDiff()
    #adtype = Optimization.AutoReverseDiff()
    #adtype = Optimization.AutoEnzyme()   #failed, don't use "not implemented yet error".

    optf = Optimization.OptimizationFunction((θ, p) -> loss(θ), adtype)

    #setup the bounds for the bathymetry
    lb_p = zeros(mesh.nCells)
    lb_p .= -0.2  

    ub_p = zeros(mesh.nCells)
    ub_p .= 0.5 

    optprob = Optimization.OptimizationProblem(optf, ps, lb=lb_p, ub=ub_p)
    #optprob = Optimization.OptimizationProblem(optf, ps)

    #res = Optimization.solve(optprob, PolyOpt(), callback = callback)  #PolyOpt does not support lb and ub 
    #res = Optimization.solve(optprob, NLopt.LD_LBFGS(), callback = callback)   #very fast 
    #res = Optimization.solve(optprob, Optim.BFGS(), callback=callback; iterations=30, maxiters=40, f_calls_limit=20, show_trace=true)
    res = Optimization.solve(optprob, Optim.BFGS(), callback=callback, maxiters = 1, local_maxiters = 3; f_tol=1e-3, show_trace=false)  #iterations=10 
    #res = Optimization.solve(optprob, Optim.LBFGS(), callback=callback)  #oscilates around 1e-7
    #res = Optimization.solve(optprob, Optim.Newton(), callback=callback)  #error: not supported as the Fminbox optimizer
    #res = Optimization.solve(optprob, Optim.GradientDescent(), callback=callback)  #very slow decrease in loss 
    #res = Optimization.solve(optprob, Flux.Adam(0.01), callback=callback, maxiters=1000)
    
    @show res
    @show res.original

    #save the inversion results
    #@show PARS 
    jldsave(joinpath(save_path, "inversion_results.jld2"); LOSS, PRED, PARS, sol)

end

if bPlot_Inversion_Results

    println("   Plotting inversion results ...")

    #plot inversion results
    plot_inversion_results(save_path, mesh, zb_cell_truth)
end 

time_end = Dates.now()

time_elapsed = Dates.canonicalize(Dates.CompoundPeriod(time_start-time_end))
println("Wall time: ", time_elapsed)

println("Done!")
