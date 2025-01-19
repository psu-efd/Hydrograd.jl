#This code solvers a 1D SWE with the following:
# BCs: inletQ on the left and exitH on the right
# Bed profile: prescribed 
# ICs: prescribed 
# 
#Then the code does an inversion of three parameters [inletQ, exitH, ManningN]

using Plots
using JLD2

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
bPerform_Inversion = false           #whether to do inversion 
bPerform_Bayesian_Estimation = true  #whether to use Bayesian inference to do the inversion (parameter estimation). 
bPlot_Inversion_Results = false     #whehter to plot the inversion results

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
    dt=0.005,          #time step size
    CFL=0.4,           #CFL number 
    tEnd=30.0,          #end of simulation time 
    h_small=1.0e-2,      #a small water depth for dry/wet treatment
    RiemannSolver="HLL"  #choice of Riemann problem solver 
)

#variables

ManningN = 0.03    #Manning's n

#  setup bed 
zb_face = zeros(mesh.nFaces)      #zb at faces (points in 1D)
zb_cell = zeros(mesh.nCells)      #zb at cell centers 
S0 = zeros(mesh.nCells)          #bed slope at cell centers 

setup_bed!(mesh, zb_face, zb_cell, S0)

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

#set up ODE parameters. In this problem, p=[inletQ, exitH, ManningN]
p = [inletQ, exitH, ManningN]
Q0 = hcat(h, q)   #initial condition: Q = [h q]

# time information 
tspan = (0.0, swe_1d_constants.tEnd)
dt = swe_1d_constants.dt
t = 0.0:dt:swe_1d_constants.tEnd

# define the ODE
f = ODEFunction((dQdt, Q, p, t) -> swe_1D_rhs!(dQdt, Q, p, t,
        mesh, swe_1d_constants, left_bcType, right_bcType, S0), jac_prototype=nothing)

prob = ODEProblem(f, Q0, tspan, p)

if bSimulate_Synthetic_Data

    println("   Performing 1D SWE simulation ...")

    #for direct sensivity analysis
    #prob = ODEForwardSensitivityProblem((Q, p, t) -> swe_1D_rhs!(Q, p, t, 
    #                          mesh, swe_1d_constants, left_bcType, right_bcType, S0), Q0, tspan, p)

    # solve the ODE with a solver of your choice 
    sol = solve(prob, Heun(), adaptive=false, dt=dt, saveat=t)
    #sol = solve(prob, AB3(), adaptive=false, dt=dt, saveat = t)
    #sol = solve(prob, Rosenbrock23(), adaptive=false, dt=dt, saveat = t)  #implicit

    #x, dp = extract_local_sensitivities(sol)
    #da = dp[1]
    #plot(sol.t, da', lw = 3)
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

    #inversion parameters: inletQ, exitH, and ManningN
    # initial guess for the parameters 
    ps = [0.5, 0.5, 0.02]

    function predict(θ)
        Array(solve(prob, Heun(), adaptive=false, p=θ, dt=dt, saveat=t)) #[:,1,end]
    end

    SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true
    #jac = ForwardDiff.jacobian(predict, ps)
    #jac = Zygote.jacobian(predict, ps)
    #jac = ReverseDiff.jacobian(predict, ps)
    #@show jac
    #plot(jac)
    #exit()

    ## Defining Loss function
    function loss(θ)
        pred = predict(θ)
        l = pred - Array(sol)
        l = l[:, 1, end]  #only take the last time step water depth result
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
        append!(PARS, [θ])

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

    adtype = Optimization.AutoZygote()
    #adtype = Optimization.AutoForwardDiff()
    #adtype = Optimization.AutoReverseDiff()
    #adtype = Optimization.AutoEnzyme()   #failed, don't use "not implemented yet error".

    optf = Optimization.OptimizationFunction((θ, p) -> loss(θ), adtype)

    optprob = Optimization.OptimizationProblem(optf, ps, lb=[0.1, 0.4, 0.01], ub=[1.0, 1.0, 0.06])

    #res = Optimization.solve(optprob, PolyOpt(), callback = callback)  #PolyOpt does not support lb and ub 
    #res = Optimization.solve(optprob, NLopt.LD_LBFGS()(), callback = callback)
    res = Optimization.solve(optprob, Optim.BFGS(), callback=callback, maxiters=100)

    @show res.u

    #save the inversion results
    jldsave(joinpath(save_path, "inversion_results.jld2"); LOSS, PRED, PARS, sol)

end

if bPerform_Bayesian_Estimation

    println("   Performing Bayesian estimation ...")

    #open the simulation result (as the ground truth)
    sol = load(joinpath(save_path, "simulation_solution.jld2"))["sol"]
    #@show sol

    #add some random noise to the observations of water depth h 
    h_with_noise = Array(sol)[:,1,end] + (rand(size(Array(sol)[:,1,end])[1]) .- 0.5)/10

    #plot simulation and noisy observations 
    p1 = plot(mesh.xCells, Array(sol)[:,1,end] + zb_cell, alpha=0.3, c=:aqua, label="simulation")
    scatter!(mesh.xCells, h_with_noise + zb_cell, mc=:blue, ms=2, ma=0.5, label="noisy data" )

    display(p1)

    @model function fitlv(data, prob)
        # Prior distributions.
        σ ~ InverseGamma(2, 3)   #how to determine the parameters?
        inletQ ~ truncated(Normal(1.0, 0.5); lower=0.1, upper=2.0)
        exitH  ~ truncated(Normal(1.0, 0.5); lower=0.1, upper=2.0)
        ManningN ~ truncated(Normal(0.02, 0.01); lower=0.01, upper=0.06)
 
    
        # Simulate SWE_1D model. 
        ps_Bayesian = [inletQ, exitH, ManningN]
        
        # Need to redefine the ODEProblem with ps_Bayesian?
        #prob_Bayesian = ODEProblem(f, Q0, tspan, ps_Bayesian)

        predicted = solve(prob, Heun(); p=ps_Bayesian, adaptive=false, dt=dt, saveat=t)
    
        # Observations.
        #for i in 1:length(predicted)
        #    data[:, i] ~ MvNormal(predicted[i], σ^2 * I)
        #end
        data[:] ~ MvNormal(predicted[:,1,end], σ^2 * I)   #only use water depth at tEnd
    
        return nothing
    end
    
    model = fitlv(h_with_noise, prob)
    
    # Sample 3 independent chains with forward-mode automatic differentiation (the default).
    chain = sample(model, NUTS(0.65), MCMCSerial(), 100, 1; progress=true)  #serial 
    #chain = sample(model, NUTS(0.65), MCMCThreads(), 100, 3; progress=true)  #multiple threads

    #save the chain 
    # Save a chain.
    write("chain-file.jls", chn)

    # Read the chain back 
    #chain_reloaded = read("chain-file.jls", chain)

    p2=plot(chain)
    display(p2)

    println("type enter to exit")
    readline()
    exit()
end 

if bPlot_Inversion_Results

    println("   Plotting inversion results ...")

    #plot inversion results
    plot_inversion_results(save_path, mesh, zb_cell)
end 

println("Done!")
