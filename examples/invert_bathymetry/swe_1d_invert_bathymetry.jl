#This code solvers a 1D SWE with the following:
# BCs: inletQ on the left and exitH on the right
# Bed profile: prescribed 
# ICs: prescribed 
# 
#Then the code does an inversion of the bathymetry

using Plots
using JLD2
using BenchmarkTools

using AdHydraulics

#SciML related packages
using DifferentialEquations, Optimization, OptimizationPolyalgorithms, OptimizationNLopt,
    Zygote, SciMLSensitivity
using ModelingToolkit
using Symbolics

using Optim

using ForwardDiff, ReverseDiff

using OrdinaryDiffEq

#for Bayesian estimation
using Turing

# Load StatsPlots for visualizations and diagnostics.
using StatsPlots

using LinearAlgebra

using Random
Random.seed!(1234)

#include problem-specific Julia files (must be in the same directory as the main file)
include("process_bed.jl")
include("process_initial_condition.jl")
include("semi_discretize.jl")
include("plot_inversion_results.jl")

#print some information
println("SWE-1D simulation and inversion ...")

#directory to save to (the same directory as the main file)
save_path = dirname(@__FILE__)

#define control variables
bSimulate_Synthetic_Data = false    #whether to do the 1D SWE simulation to create synthetic data 
bPlot_Simulation_Results = false    #whether to plot simulation results
bPerform_Inversion = true           #whether to do inversion 
bPlot_Inversion_Results = true     #whehter to plot the inversion results

#mesh related variables
nCells = 100     # Number of cells
L = 10.0  # Length of the domain

#create the 1D mesh
mesh = initialize_mesh_1D(nCells, L)

#create problem constants 
swe_1d_constants = swe_1D_const(
    g=9.81,
    t=0.0,             #starting time (current time)
    dt_min=0.005,      #minimum dt (if adaptive time step)
    dt=0.001,          #time step size
    CFL=0.4,           #CFL number 
    tEnd=20.0,          #end of simulation time 
    h_small=1.0e-3,      #a small water depth for dry/wet treatment
    RiemannSolver="HLL"  #choice of Riemann problem solver 
)

#variables

ManningN = 0.03    #Manning's n

#  setup bed 
zb_face = zeros(Number, mesh.nFaces)      #zb at faces (points in 1D)
zb_cell = zeros(Number, mesh.nCells)      #zb at cell centers 
S0 = zeros(Number, mesh.nCells)          #bed slope at cell centers 

setup_bed!(mesh, zb_face, zb_cell, S0)

#make a copy of the bathymetry truth: zb_face_truth, zb_cell_truth
zb_face_truth = deepcopy(zb_face)
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
        mesh, swe_1d_constants, left_bcType, right_bcType, left_bcValue, right_bcValue, ManningN,
        zb_face), 
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
    swe_1D_make_animation(sol, mesh, zb_cell, save_path)
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
        Array(solve(prob, Tsit5(), adaptive=true, dt=dt, saveat=t))[:,1,end]
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
        l = pred - Array(sol)[:,1,end] + (θ - zb_cell_truth)  #loss = free surface elevation mismatch

        return sum(abs2, l), pred # Mean squared error
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

    callback = function (θ, l, pred) #callback function to observe training
        iter = size(LOSS)[1]  #get the inversion iteration number (=length of LOSS array)
        println("      iter, loss = ", iter, ", ", l)

        append!(PRED, [pred])
        append!(LOSS, l)

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

    #adtype = Optimization.AutoZygote()         #works, but slow. Maybe because of memory allocations?
    adtype = Optimization.AutoForwardDiff()
    #adtype = Optimization.AutoReverseDiff()
    #adtype = Optimization.AutoEnzyme()   #failed, don't use "not implemented yet error".

    optf = Optimization.OptimizationFunction((θ, p) -> loss(θ), adtype)

    #setup the bounds for the bathymetry
    lb_p = zeros(mesh.nCells)
    lb_p .= -0.2  

    ub_p = zeros(mesh.nCells)
    ub_p .= 0.5 

    optprob = Optimization.OptimizationProblem(optf, ps, lb=lb_p, ub=ub_p)

    #res = Optimization.solve(optprob, PolyOpt(), callback = callback)  #PolyOpt does not support lb and ub 
    #res = Optimization.solve(optprob, NLopt.LD_LBFGS(), callback = callback)   #very fast 
    #res = Optimization.solve(optprob, Optim.BFGS(), callback=callback; iterations=30, maxiters=40, f_calls_limit=20, show_trace=true)
    res = Optimization.solve(optprob, Optim.BFGS(), callback=callback)
    #res = Optimization.solve(optprob, Optim.LBFGS(), callback=callback)  #oscilates around 1e-7
    #res = Optimization.solve(optprob, Optim.Newton(), callback=callback)  #error: not supported as the Fminbox optimizer
    #res = Optimization.solve(optprob, Optim.GradientDescent(), callback=callback)  #very slow decrease in loss 
    
    @show res.u

    #save the inversion results
    #@show PARS 
    jldsave(joinpath(save_path, "inversion_results.jld2"); LOSS, PRED, PARS, sol)

end

if bPlot_Inversion_Results

    println("   Plotting inversion results ...")

    #plot inversion results
    plot_inversion_results(save_path, mesh, zb_face, zb_cell, S0, zb_cell_truth)
end 

println("Done!")
