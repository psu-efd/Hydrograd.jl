#This code solvers a 1D SWE with the following:
# BCs: inletQ on the left and exitH on the right
# Bed profile: prescribed 
# ICs: prescribed 
# 
#Then the code does an inversion of three parameters [inletQ, exitH, ManningN]

using Plots

using AdHydraulics

#SciML related packages
using DifferentialEquations, Optimization, OptimizationPolyalgorithms, OptimizationNLopt,
    Zygote, SciMLSensitivity
using ModelingToolkit
using Symbolics

using Optim 

using ForwardDiff, ReverseDiff

using OrdinaryDiffEq


include("process_bed.jl")
include("process_initial_condition.jl")
include("semi_discretize.jl")


#define variables

#mesh related 
nCells = 20     # Number of cells
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
    tEnd=0.1,          #end of simulation time 
    h_small=1.0e-3,      #a small water depth for dry/wet treatment
    RiemannSolver="HLL"  #choice of Riemann problem solver 
)

#variables

ManningN = 0.03    #Manning's n

#  setup bed 
zb_face = zeros(mesh.nFaces)      #zb at faces (points in 1D)
zb_cell = zeros(mesh.nCells)      #zb at cell centers 
S0 =  zeros(mesh.nCells)          #bed slope at cell centers 

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
left_bcValue = inletQ = 0.5
right_bcType = "exitH"
#right_bcType = "wall"
right_bcValue = exitH = 1.0

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

#for direct sensivity analysis
#prob = ODEForwardSensitivityProblem((Q, p, t) -> swe_1D_rhs!(Q, p, t, 
#                          mesh, swe_1d_constants, left_bcType, right_bcType, S0), Q0, tspan, p)

# solve the ODE with a solver of your choice 
sol = solve(prob, Heun(), adaptive=false, dt = dt, saveat = t) 
#sol = solve(prob, AB3(), adaptive=false, dt=dt, saveat = t)
#sol = solve(prob, Rosenbrock23(), adaptive=false, dt=dt, saveat = t)  #implicit

#x, dp = extract_local_sensitivities(sol)
#da = dp[1]
#plot(sol.t, da', lw = 3)
#readline()
#exit()

swe_1D_save_results(sol, total_water_volume, mesh, zb_cell)
swe_1D_make_plots()
#readline()
#exit()

################
#Inversion part# 
################

println("Entering inversion ...")

#inverse the parameters: inletQ, exitH, and ManningN
# initial guess for the parameters 
ps = [0.2, 0.8, 0.02]

function predict(θ)
    Array(solve(prob, Heun(), adaptive=false, p = θ, dt = dt, saveat = t)) #[:,1,end]
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
    l = l[:,1,end]  #only take the last time step water depth result
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
    println("loss = ", l)
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

optprob = Optimization.OptimizationProblem(optf, ps, lb = [0.1, 0.4, 0.01], ub = [1.0, 1.0, 0.06])

#res = Optimization.solve(optprob, PolyOpt(), callback = callback)  #PolyOpt does not support lb and ub 
#res = Optimization.solve(optprob, NLopt.LD_LBFGS()(), callback = callback)
res = Optimization.solve(optprob, Optim.BFGS(), callback = callback)

@show res.u 

println("Done!")
