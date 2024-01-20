using Plots
using CSV
using DataFrames
using Printf

#using AdHydraulics

#SciML related packages
using DifferentialEquations, Optimization, OptimizationPolyalgorithms, OptimizationNLopt,
    Zygote, SciMLSensitivity
using ModelingToolkit
using Symbolics

using Optim 

using ForwardDiff, ReverseDiff

using BandedMatrices

using OrdinaryDiffEq

#SWE 1D solver all in one file and simplified

#Functions to setup the bed profile (this is just for some examples)
function my_gauss(x::Float64; sigma::Float64=1.0, h::Float64=1.0, mid::Float64=0.0)
    variance = sigma^2
    return h * exp(-(x - mid)^2 / (2 * variance))
end

function setup_bed!(nCells, nFaces, L, xFaces, zb_face, zb_cell, S0)
    # parameters for bed setup
    b_bump_height = 0.3

    #loop through faces (points)
    @inbounds for i in 1:nFaces
        zb_face[i] = my_gauss(xFaces[i], sigma=1.0, h=b_bump_height, mid=L / 2.0)
    end

    interploate_zb_from_face_to_cell_and_compute_S0!(nCells, xFaces, zb_face, zb_cell, S0)

    # optionally plot the bed for checking
    plot_bed(xCells, zb_cell, zb_face, S0)
end

function interploate_zb_from_face_to_cell_and_compute_S0!(nCells, xFaces, zb_face, zb_cell, S0)
    #loop through cells
    @inbounds for i in 1:nCells
        zb_cell[i] = (zb_face[i+1] + zb_face[i]) / 2.0
        S0[i] = -(zb_face[i+1] - zb_face[i]) / (xFaces[i+1] - xFaces[i])
    end
end

#plot the bed profile for checking 
function plot_bed(xCells, zb_cell, zb_face, S0)
    p1 = plot(xCells, zb_cell, label="bed at cells")
    plot!(xFaces, zb_face, label="bed at faces")
    p2 = plot(xCells, S0, label="slope at cells")
    display(p1)
    display(p2)
end

# setup initial condition: free surface 
function setup_initial_eta!(nCells, xCells, eta, zb_cell, h, h_small)
    # parameters for initial free surface profile setup
    bump_center_x = 5.0  # center of the bump
    bump_half_width_x = 1.0  # bump's half width

    #loop over cells
    @inbounds for i in 1:nCells
        #if abs(mesh.xCells[i] - bump_center_x) < bump_half_width_x:
        #    fields.eta[i] = 1.1
        #else:
        #    fields.eta[i] = 1.0

        if xCells[i] < bump_center_x
            eta[i] = 1.0
        else
            eta[i] = 0.5 #0.5  #0.5
        end
    end

    #update water depth
    @inbounds for i in 1:nCells
        h[i] = eta[i] - zb_cell[i]
    end

    h[h.<0.0] .= h_small  #ensure positivity of water depth h

    #update the free surface elevation again in case h has been clipped
    @inbounds for i in 1:nCells
        eta[i] = h[i] + zb_cell[i]
    end

    #optionally plot the free surface for checking 
    plot_free_surface_elevation(xCells, zb_cell, eta, h)
end

#plot the free surface profile for checking 
function plot_free_surface_elevation(xCells, zb_cell, eta, h)
    p1 = plot(xCells, zb_cell, label="bed at cells")
    plot!(xCells, eta, label="free surface")
    plot!(xCells, h, label="water depth")
    display(p1)
end

# HLL Riemman solver
#Ref: https://github.com/wme7/ShallowWaterEquations
# bcType: -1 (for internal faces, not even a boundary face)
function Riemann_1D_hll(g, h_face, q_face, bcType, bcValue, h_small)
    h_L = h_face[1]
    h_R = h_face[2]
    q_L = q_face[1]
    q_R = q_face[2]

    if (h_L <= h_small && h_R <= h_small)
        if bcType=="inletQ"    #for inletQ BC, we need to specify the discharge (q=hu) even it is dry
            return [bcValue, zero(eltype(bcValue))]
        else
            return [zero(eltype(bcValue)), zero(eltype(bcValue))]
        end 
    end 

    if h_L <= h_small
        u_L = zero(eltype(bcValue))
    else
        u_L = q_L / h_L
    end

    if h_R <= h_small
        u_R = zero(eltype(bcValue))
    else
        u_R = q_R / h_R
    end

    if h_L <= h_small
        s_L = u_R - 2.0 * sqrt(g * h_R)
    else
        s_L = min(u_R - sqrt(g * h_R), u_L - sqrt(g * h_L))
    end 

    if h_R <= h_small
        s_R = u_L + 2.0 * sqrt(g * h_L)
    else
        s_R = max(u_R + sqrt(g * h_R), u_L + sqrt(g * h_L))
    end

    h_flux_L = q_L
    h_flux_R = q_R

    if (abs(s_R-s_L) < 1e-10)
        h_flux_star = zero(eltype(bcValue))
    else
        h_flux_star = (s_R * h_flux_L - s_L * h_flux_R + s_L * s_R * (h_R - h_L)) / (s_R - s_L)
    end 

    q_flux_L = q_L * u_L + g * h_L * h_L / 2.0
    q_flux_R = q_R * u_R + g * h_R * h_R / 2.0

    if (abs(s_R-s_L) < 1e-10)
        q_flux_star = zero(eltype(bcValue))
    else
        q_flux_star = (s_R * q_flux_L - s_L * q_flux_R + s_L * s_R * (q_R - q_L)) / (s_R - s_L)
    end 

    if (0.0 <= s_L)
        h_flux = h_flux_L
        q_flux = q_flux_L
    elseif (s_L <= 0 && 0 <= s_R)
        h_flux = h_flux_star
        q_flux = q_flux_star
    else
        h_flux = h_flux_R
        q_flux = q_flux_R
    end 

    return [h_flux, q_flux]
end 

# function to calculate the RHS of ODE: dQdt (= dhdt and dqdt)
function swe_1D_rhs!(dQdt, Q, p, t, dx, g, h_small, RiemannSolver,
    nCells, nFaces, left_bcType, right_bcType, S0)

    #fluxex on faces 
    fluxes = zeros(eltype(Q), 2, nFaces)

    #println("time t = ", t)

    #set the parameter values 
    left_bcValue = p[1] #hard-coded here. Should be more elegant. 
    right_bcValue = p[2]
    ManningN = p[3]

    #get a view (no need for a copy) of the current solution variables h and q 
    h = Q[:,1]
    q = Q[:,2]

    #ghost cells on the left and right boundaries
    #    wall (no flux, dh/dx=0 and q=hu=0), bcValue = 0.0 (not used)
    #    zeroGradient (dh/dx=0, dq/dx=0), bcValue = 0.0 (not used)
    #    inletQ (specify inlet specific discharge q=hu=q_in), bcValue = q_in
    #    exitH (specify outlet water depth, h = h_outlet), bcValue = h_outlet

    #left boundary
    if left_bcType=="wall"   #wall
        h_ghostCell_left = h[1]
        q_ghostCell_left = -q[1]    #flux goes the opposite direction -> no flux
    elseif left_bcType=="zeroGradient #zeroGradient"
        h_ghostCell_left = h[1]
        q_ghostCell_left = q[1]
    elseif left_bcType=="inletQ" #inletQ (bcValue=q_in)
        h_ghostCell_left = h[1]
        q_ghostCell_left = left_bcValue
    elseif left_bcType=="exitH" #exitH (bcValue=h_outlet)
        h_ghostCell_left = left_bcValue
        q_ghostCell_left = q[1]
    else 
        println("Left boundary condition type not recognized.")
        exit(-1)  #exit with an error code of -1
    end

    # right boundary
    if right_bcType == "wall"  # wall
        h_ghostCell_right = h[nCells]
        q_ghostCell_right = -q[nCells]
    elseif right_bcType == "zeroGradient"  # zeroGradient
        h_ghostCell_right = h[nCells]
        q_ghostCell_right = q[nCells]
    elseif right_bcType == "inletQ"  # inletQ (bcValue=q_in)
        h_ghostCell_right = h[nCells]
        q_ghostCell_right = right_bcValue
    elseif right_bcType == "exitH"  # exitH (bcValue=h_outlet)
        h_ghostCell_right = right_bcValue
        q_ghostCell_right = q[nCells]
    else 
        println("Right boundary condition type not recognized.")
        exit(-1)  #exit with an error code of -1
    end

    #loop through all faces
    @inbounds for iFace in 1:nFaces
        if (iFace == 1)               #left boundary face
            h_face = [h_ghostCell_left, h[1]]
            q_face = [q_ghostCell_left, q[1]]
            bcType = left_bcType
            bcValue = left_bcValue
        elseif (iFace == nFaces)   #right boundary face
            h_face = [h[nCells], h_ghostCell_right]
            q_face = [q[nCells], q_ghostCell_right]
            bcType = left_bcType
            bcValue = left_bcValue
        else                          #internal face
            h_face = [h[iFace-1], h[iFace]]
            q_face = [q[iFace-1], q[iFace]]
            bcType = "internal"
            bcValue = zero(eltype(p))  #not used
        end 

        if (RiemannSolver == "Roe")
            println("Not implemented yet")
            exit(-1)  #exit with an error code of -1
        elseif (RiemannSolver == "HLL")
            flux = Riemann_1D_hll(g, h_face, q_face, bcType, bcValue, h_small)
        elseif (RiemannSolver == "HLLC")
            println("Not implemented yet")
            exit(-1)  #exit with an error code of -1
        else
            println("Wrong choice of RiemannSolver")
            exit(-1)  #exit with an error code of -1
        end

        fluxes[1,iFace] = flux[1] 
        fluxes[2,iFace] = flux[2]
    end 

    #loop through all cells
    @inbounds for iCell in 1:nCells

        #calcuate the RHS of the ODEs
        #flux_east = @view fluxes[:,iCell+1] #flux on east face
        #flux_west = @view fluxes[:,iCell]   #flux on west face
        flux_east = fluxes[:,iCell+1] #flux on east face
        flux_west = fluxes[:,iCell]   #flux on west face

        #fields.dhdt[iCell] = - (flux_east[1]- flux_west[1]) / dx
        dQdt[iCell,1] = - (flux_east[1]- flux_west[1]) / dx

        if (h[iCell] <= h_small) #if a dry cell, no flow resistance term
            #fields.dqdt[iCell] = -(flux_east[2] - flux_west[2]) / dx + g * h[iCell] * fields.S0[iCell]
            dQdt[iCell,2] = -(flux_east[2] - flux_west[2]) / dx + g * h[iCell] * S0[iCell]
        else
            #fields.dqdt[iCell] = (- (flux_east[2] - flux_west[2]) / dx + g * h[iCell] * fields.S0[iCell]
            #           - g*ManningN^2/max(h[iCell], h_small)^(7.0/3.0)*abs(fields.q[iCell])*fields.q[iCell])
            dQdt[iCell,2] = (- (flux_east[2] - flux_west[2]) / dx + g * h[iCell] * S0[iCell]
                       - g*ManningN^2/max(h[iCell], h_small)^(7.0/3.0)*abs(q[iCell])*q[iCell])
        end
    end 

end

#calculate the total water volume in the domain
function calc_total_water_volume(h, dx)
    total_water_volume = sum(h) * dx
    return total_water_volume
end

#save results 
function save_results(sol, total_water_volume, dx, nCells, zb_cell)
    #solution at the end of the simulation
    sol_final = sol.u[end]

    #calculate total volume of water 
    for (Q, t) in tuples(sol)
        push!(total_water_volume, calc_total_water_volume(Q[:, 1], dx))
    end

    #save results to files the same directory of this jl file
    open(joinpath(dirname(@__FILE__), "results.csv"), "w") do fo
        println(fo, "x, h, q, eta, zb")
        for i in 1:nCells
            @printf(fo, "%.6f, %.6f, %.6f, %.6f, %.6f\n", xCells[i], sol_final[i, 1], sol_final[i, 2], sol_final[i, 1] + zb_cell[i], zb_cell[i])
        end
    end

    open(joinpath(dirname(@__FILE__), "total_water_volume.csv"), "w") do fo
        println(fo, "total_water_volume")
        for volume in total_water_volume
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

#define some Variables

#mesh related 

nCells = 50     # Number of cells
nFaces = nCells + 1     # Number of faces (= nCells + 1)
L = 10.0  # Length of the domain

#cell size
dx = L / nCells

#x coordinates of cell centers 
xCells = [ dx/2.0 + dx*i for i=0:(nCells-1) ]

#x coordinates of faces (points in 1D) 
xFaces = [ dx*i for i=0:(nFaces-1) ]


#solution parameters
g = 9.81   #gravity constant 

ManningN = 0.03    #Manning's n

t = 0.0   #time 
dt_min = 0.005  #minimum time step size 
dt = 0.005   #time step size 
CFL = 0.4  #CFL number 
tEnd = 1.0  #end time for simulation 

h_small = 1.0e-3  #dry bed water depth threshold (e.g., 1.0e-3)

RiemannSolver = "HLL"  #The choice of Riemann solver, e.g., Roe, HLL, and HLLC 

#variables

zb_face = zeros(nFaces)
zb_cell = zeros(nCells)
S0 =  zeros(nCells)

setup_bed!(nCells, nFaces, L, xFaces, zb_face, zb_cell, S0)

eta = zeros(nCells)
h = zeros(nCells)
q = zeros(nCells)
total_water_volume = []

setup_initial_eta!(nCells, xCells, eta, zb_cell, h, h_small)

#boundary conditions
left_bcType = "inletQ"
#left_bcType = "wall"
left_bcValue =  0.5
right_bcType = "exitH"
#right_bcType = "wall"
right_bcValue = 1.0


#set up ODE 
p = [0.5, 1.0, 0.03]
Q0 = hcat(h, q)   #initial condition: Q = [h q]

tspan = (0.0, tEnd)
dt = 0.005
t = 0.0:dt:tEnd

#sparsity = BandedMatrix(Ones(nCells,nCells), (1,1))
#f = ODEFunction((dQdt, Q, p, t) -> swe_1D_rhs!(dQdt, Q, p, t, dx, g, h_small, RiemannSolver,
#nCells, nFaces, left_bcType, right_bcType, S0), jac_prototype=sparsity)

f = ODEFunction((dQdt, Q, p, t) -> swe_1D_rhs!(dQdt, Q, p, t, dx, g, h_small, RiemannSolver,
nCells, nFaces, left_bcType, right_bcType, S0), jac_prototype=nothing)

prob = ODEProblem(f, Q0, tspan, p)

#try to use Symbolics to automatically determine Jacobian sparsity (failed)
#dQ0 = copy(Q0)
#jac_sparsity = Symbolics.jacobian_sparsity((dQdt, Q) -> swe_1D_rhs!(dQdt, Q, p, t, dx, g, h_small, RiemannSolver,
#           nCells, nFaces, left_bcType, right_bcType, S0),
#           dQ0, Q0)
#@show jac_sparsity

#sys = modelingtoolkitize(prob);
#sparseprob = ODEProblem(sys, Pair[], t, jac = true, sparse = true)

#for direct sensivity analysis
#prob = ODEForwardSensitivityProblem((Q, p, t) -> swe_1D_rhs!(Q, p, t, dx, g, h_small, RiemannSolver,
#       nCells, nFaces, left_bcType, right_bcType, S0), Q0, tspan, p)

sol = solve(prob, Heun(), adaptive=false, dt = dt, saveat = t) 
#sol = solve(prob, AB3(), adaptive=false, dt=dt, saveat = t)
#sol = solve(prob, Rosenbrock23(), adaptive=false, dt=dt, saveat = t)  #implicit

#x, dp = extract_local_sensitivities(sol)
#da = dp[1]
#plot(sol.t, da', lw = 3)
#readline()
#exit()

save_results(sol, total_water_volume, dx, nCells, zb_cell)
make_plots()
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

#l, pred = loss(ps)
#println("l, typeof(l) = ", l, " ", typeof(l))

#grad = Zygote.gradient(loss, ps)
#grad = ForwardDiff.gradient(loss, ps)
#grad = ForwardDiff.gradient(predict, θ)
#@show grad
#println("grad = ", grad)
#readline()
#exit()


#l, pred = loss(ps)
#size(pred), size(sol), size(t) # Checking sizes

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
