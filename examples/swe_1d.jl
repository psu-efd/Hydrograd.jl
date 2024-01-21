"""
This is a simple example for solving 1D SWE using Godnove scheme. 

The development of this 1D referenced the following Python implementation:
Ref: https://github.com/wme7/ShallowWaterEquations
"""

using Plots
using CSV
using DataFrames
using Printf

using AdHydraulics

#SciML related packages
using DifferentialEquations, Optimization, OptimizationPolyalgorithms, OptimizationNLopt,
    Zygote, SciMLSensitivity
using ModelingToolkit
using Symbolics
using Enzyme

using Optim 

using ForwardDiff, ReverseDiff

using BandedMatrices

using OrdinaryDiffEq



#Functions to setup the bed profile (this is just for some examples)
function my_gauss(x::Float64; sigma::Float64=1.0, h::Float64=1.0, mid::Float64=0.0)
    variance = sigma^2
    return h * exp(-(x - mid)^2 / (2 * variance))
end

function setup_bed!(mesh::mesh_1D, fields::swe_1D_fields)
    # parameters for bed setup
    b_bump_height = 0.3

    #loop through faces (points)
    for i in 1:mesh.nFaces
        fields.zb_face[i] = my_gauss(mesh.xFaces[i], sigma=1.0, h=b_bump_height, mid=mesh.L / 2.0)
    end

    #loop through cells
    for i in 1:mesh.nCells
        fields.zb_cell[i] = (fields.zb_face[i+1] + fields.zb_face[i]) / 2.0
        fields.S0[i] = -(fields.zb_face[i+1] - fields.zb_face[i]) / (mesh.xFaces[i+1] - mesh.xFaces[i])
    end

    # optionally plot the bed for checking
    #plot_bed(mesh, fields)
end

#plot the bed profile for checking 
function plot_bed(mesh::mesh_1D, fields::swe_1D_fields)
    p1 = plot(mesh.xCells, fields.zb_cell, label="bed at cells")
    plot!(mesh.xFaces, fields.zb_face, label="bed at faces")
    p2 = plot(mesh.xCells, fields.S0, label="slope at cells")
    display(p1)
    display(p2)
end

# setup initial condition: free surface 
function setup_initial_eta!(mesh, fields)
    # parameters for initial free surface profile setup
    bump_center_x = 5.0  # center of the bump
    bump_half_width_x = 1.0  # bump's half width

    #loop over cells
    for i in 1:mesh.nCells
        #if abs(mesh.xCells[i] - bump_center_x) < bump_half_width_x:
        #    fields.eta[i] = 1.1
        #else:
        #    fields.eta[i] = 1.0

        if mesh.xCells[i] < bump_center_x
            fields.eta[i] = 1.0
        else
            fields.eta[i] = 0.5 #0.0  #0.5
        end
    end

    #update water depth
    fields.h = fields.eta - fields.zb_cell
    fields.h[fields.h.<0.0] .= fields.swe_1d_para.h_small  #ensure positivity of water depth h

    #update the free surface elevation again in case h has been clipped
    fields.eta = fields.h + fields.zb_cell

    #optionally plot the free surface for checking 
    #plot_free_surface_elevation(mesh, fields)
end

#plot the free surface profile for checking 
function plot_free_surface_elevation()
    p1 = plot(mesh.xCells, fields.zb_cell, label="bed at cells")
    plot!(mesh.xCells, fields.eta, label="free surface")
    plot!(mesh.xCells, fields.h, label="water depth")
    display(p1)
end

#calculate the total water volume in the domain
function calc_total_water_volume(h, dx)
    total_water_volume = sum(h) * dx
    return total_water_volume
end

#save results 
function save_results(sol, mesh, fields)
    #solution at the end of the simulation
    sol_final = sol.u[end]

    #calculate total volume of water 
    for (Q, t) in tuples(sol)
        push!(fields.total_water_volume, calc_total_water_volume(Q[:, 1], mesh.dx))
    end

    #save results to files the same directory of this jl file
    open(joinpath(dirname(@__FILE__), "results.csv"), "w") do fo
        println(fo, "x, h, q, eta, zb")
        for i in 1:mesh.nCells
            @printf(fo, "%.6f, %.6f, %.6f, %.6f, %.6f\n", mesh.xCells[i], sol_final[i, 1], sol_final[i, 2], sol_final[i, 1] + fields.zb_cell[i], fields.zb_cell[i])
        end
    end

    open(joinpath(dirname(@__FILE__), "total_water_volume.csv"), "w") do fo
        println(fo, "total_water_volume")
        for volume in fields.total_water_volume
            println(fo, volume)
        end
    end

end

#make plots of the results
function make_plots()
    data_dataframe = CSV.read(joinpath(dirname(@__FILE__), "results.csv"), DataFrame)
    data = Matrix(data_dataframe)  #x, h, q, eta, zb

    p1 = plot(data[:, 1], data[:, 4], label="free surface")
    plot!(data[:, 1], data[:, 5], label="bed")
    plot!(legend=:outerbottom, fg_legend=:transparent)

    data_dataframe = CSV.read(joinpath(dirname(@__FILE__), "total_water_volume.csv"), DataFrame)
    data = Matrix(data_dataframe)  #total_water_volume

    p2 = plot(data[:, 1])

    display(p1)
    #display(p2)
end

#make animation
function make_animation(sol, mesh, fields)

    anim = Animation()

    let animCounter = 1

        for (Q, t) in tuples(sol)

            animCounter = animCounter + 1

            if (animCounter % 10) == 0

                plt = plot(mesh.xCells, Q[:, 1] + fields.zb_cell, label="free surface")
                plot!(mesh.xCells, fields.zb_cell, label="bed")
                xlims!(-0.1, 10.1)
                ylims!(-0.1, 2.0)
                plot!(legend=:bottomleft, fg_legend=:transparent)
                frame(anim, plt)

            end

        end
    end

    #gif(anim, joinpath(dirname(@__FILE__), "free_surface_1d.gif"), fps=15)
    gif(anim, joinpath(dirname(@__FILE__), "free_surface_1d.mp4"), fps=15)
end

println("Starting SWE_1D ...")

#create the 1D FVM mesh 
L = 10.0 #length of the domain
nCells = 50 #number of cells

mesh = initialize_mesh_1D(nCells, L)

#println(mesh.xCells)

#create 1D FVM fields 
#  create parameters
swe_1d_para = swe_1D_const(
    g=9.81,
    t=0.0,
    dt_min=0.005,
    dt=0.005,
    CFL=0.4,
    tEnd=1.0,
    h_small=1.0e-3,
    RiemannSolver="HLL"
)

#  create the left and right boundaries
#boundary conditions
left_bcType = "inletQ"
#left_bcType = "wall"
left_bcValue =  0.5
right_bcType = "exitH"
#right_bcType = "wall"
right_bcValue = 1.0

fields = initialize_swe_1D_fields(mesh)

# setup the bed 
setup_bed!(mesh, fields)

#setup the initial free surface 
setup_initial_eta!(mesh, fields)

ManningN=0.03,

#Solution parameters p: inletQ, exitH, and ManningN 
p = [left_bcValue, right_bcValue, ManningN]
Q0 = hcat(fields.h, fields.q)   #initial condition: Q = [h q]

tspan = (0.0, fields.swe_1d_para.tEnd)
dt = 0.005
tEnd = fields.swe_1d_para.tEnd
t = 0.0:dt:tEnd

# solve on PDE with the Method of Line
# define the ODEProblem
#prob = ODEProblem(swe_1D_rhs!, Q0, tspan, p)
#swe_1D_rhs_with_extra(Q, p, t) = swe_1D_rhs!(Q, p, t, mesh, fields)
#prob = ODEProblem((Q, p, t) -> swe_1D_rhs!(Q, p, t, mesh, fields), Q0, tspan, p)

f = ODEFunction((dQdt, Q, p, t) -> swe_1D_rhs!(dQdt, Q, p, t, mesh, fields), jac_prototype=nothing) #closure
prob = ODEProblem(f, Q0, tspan, p)

# solve with a ODE solver of choice 
sol = solve(prob, Heun(), adaptive=false, dt = dt, saveat = t) 
#sol = solve(prob, Heun(), dt = dt, saveat = t) 
#sol = solve(prob, Tsit5(), dt = dt, saveat = t)  #Tsit5 diverged
#sol = solve(prob, Rosenbrock23(autodiff=false), dt = dt, saveat = t) #with autodiff->diverged. without AD, time step is very small.
#sol = solve(prob, DP5(), dt = dt, saveat = t) 
#sol = solve(prob, SSPRK22(), dt = dt, saveat = t) 
#sol = solve(prob, SSPRK33(), dt=dt, saveat=tEnd)

#save the results: calculate the total water volume, save final result to files 
save_results(sol, mesh, fields)

#make plots
make_plots()

#make animation 
#make_animation(sol, mesh, fields)

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
    Array(solve(prob, Heun(), adaptive=false, p = θ, dt = dt, saveat = t))
end

SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true

## Defining Loss function
function loss(θ)
    pred = predict(θ)
    l = pred - Array(sol)
    l = l[:,1,end]  #only take the last time step water depth result
    return sum(abs2, l), pred # Mean squared error
end

#l, pred = loss(ps)
#println("l, typeof(l) = ", l, " ", typeof(l))

#grad = Zygote.gradient(loss, ps)
#grad = ForwardDiff.gradient(loss, ps)

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
    #if l > 1e-3
    #    false
    #else 
    #    true
    #end

    false
end

#callback(ps, loss(ps)...) # Testing callback function

# Let see prediction vs. Truth
#scatter(sol[:, end], label = "Truth", size = (800, 500))
#plot!(PRED[end][:, end], lw = 2, label = "Prediction")

#AutoForwardDiff(): The fastest choice for small optimizations
#AutoReverseDiff(compile=false): A fast choice for large scalar optimizations
#AutoTracker(): Like ReverseDiff but GPU-compatible
#AutoZygote(): The fastest choice for non-mutating array-based (BLAS) functions
#AutoFiniteDiff(): Finite differencing, not optimal but always applicable
#AutoModelingToolkit(): The fastest choice for large scalar optimizations
#AutoEnzyme(): Highly performant AD choice for type stable and optimized code

#Enzyme.API.runtimeActivity!(true) #to get rid of the "Enzyme execution failed warning
adtype = Optimization.AutoZygote()
#adtype = Optimization.AutoForwardDiff()
#adtype = Optimization.AutoReverseDiff()
#adtype = Optimization.AutoEnzyme()  #failed, don't use for now. 

optf = Optimization.OptimizationFunction((θ, p) -> loss(θ), adtype)

optprob = Optimization.OptimizationProblem(optf, ps, lb = [0.1, 0.4, 0.01], ub = [1.0, 1.0, 0.06])

#res = Optimization.solve(optprob, PolyOpt(), callback = callback)  #PolyOpt does not support lb and ub 
#res = Optimization.solve(optprob, NLopt.LD_LBFGS()(), callback = callback)
res = Optimization.solve(optprob, Optim.BFGS(), callback = callback)

@show res.u 

println("Done!")
