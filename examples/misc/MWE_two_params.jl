using SciMLSensitivity
using DelimitedFiles, Plots
using OrdinaryDiffEq, Optimization
using OptimizationOptimisers, OptimizationPolyalgorithms, Zygote

# Problem setup parameters:
Lx = 10.0
x = 0.0:0.1:Lx
dx = x[2] - x[1]
Nx = size(x)
@show Nx

u0 = exp.(-(x .- 3.0) .^ 2) # I.C
@show u0

## Problem Parameters
p = [1.0, 1.0, 0.5]    # True solution parameters
xtrs = [dx, Nx]      # Extra parameters
dt = 0.1 * dx^2    # CFL condition
t0, tMax = 0.0, 10 * dt
tspan = (t0, tMax)
t = t0:dt:tMax;

#create material ID for each grid point
matIDs = [i <= Nx[1] ÷ 2 ? 1 : 2 for i in 1:Nx[1]]
@show matIDs

## Definition of Auxiliary functions
function ddx(u, dx)
    """
    2nd order Central difference for 1st degree derivative
    """
    return [[zero(eltype(u))]; (u[3:end] - u[1:(end - 2)]) ./ (2.0 * dx); [zero(eltype(u))]]
end

function d2dx(u, dx)
    """
    2nd order Central difference for 2nd degree derivative
    """
    return [zero(eltype(u));
            (@view(u[3:end]) .- 2.0 .* @view(u[2:(end - 1)]) .+ @view(u[1:(end - 2)])) ./
            (dx^2)
            zero(eltype(u))]
end

function update_a12(p, xtrs)
    # Model parameters
    dx, Nx = xtrs #[1.0,3.0,0.125,100]

    #the diffusion coefficient is a1 on the left half and a2 on the right half
    a12 = [p[matIDs[i]+1] for i in 1:Nx[1]]
    return a12
end

## ODE description of the Physics:
function heat(u, p, t, xtrs)
    # Model parameters
    #a0, a1, a2 = p
    dx, Nx = xtrs #[1.0,3.0,0.125,100]

    #the diffusion coefficient is a1 on the left half and a2 on the right half
    a12 = update_a12(p, xtrs)

    Zygote.ignore() do
        #println("a12: ", a12)
        #println(" ")
    end

    rhs = map(1:Nx[1]) do i  
        2.0 * p[1] * u[i] + a12[i]^2 * d2dx(u, dx)[i]
    end

    return rhs    
end
heat_closure(u, p, t) = heat(u, p, t, xtrs)

# Testing Solver on linear PDE
prob = ODEProblem(heat_closure, u0, tspan, p)
sol = solve(prob, Tsit5(), dt = dt, saveat = t);
arr_sol = Array(sol)

@show size(arr_sol)
@show arr_sol[:,end]

ps = [0.35, 0.512, 0.328];   # Initial guess for model parameters
function predict(θ)
    Array(solve(prob, Tsit5(), p = θ, dt = dt, saveat = t))[:,end]   #only take the last time step
end

## Defining Loss function
function loss(θ)
    pred = predict(θ)
    return sum(abs2.(pred .- arr_sol[:,end])) # Mean squared error
end

LOSS = []                              # Loss accumulator
PRED = []                              # prediction accumulator
PARS = []                              # parameters accumulator

cb = function (st, l) #callback function to observe training
    Zygote.ignore() do
        iter_num = length(LOSS) + 1
        println("iter_num: ", iter_num, " loss: ", l, " parameters: ", st.u)
        append!(LOSS, l)
        append!(PARS, [st.u])
    end
    
    false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)

optprob = Optimization.OptimizationProblem(optf, ps)

#res = Optimization.solve(optprob, PolyOpt(), callback = cb, maxiters = 10)
res = solve(optprob, Adam(0.01), maxiters=2, callback = cb)

@show res.u 