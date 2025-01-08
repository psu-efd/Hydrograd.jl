# UDE Example: https://docs.sciml.ai/Overview/stable/showcase/missing_physics/

# SciML Tools
using OrdinaryDiffEq, ModelingToolkit, DataDrivenDiffEq, SciMLSensitivity, DataDrivenSparse
using Optimization, OptimizationOptimisers, OptimizationOptimJL, LineSearches

using ComponentArrays   

# Standard Libraries
using LinearAlgebra, Statistics

# External Libraries
using Lux, Zygote, StableRNGs

using Plots
gr() #using GR backend
#or
#using PyPlot
#pyplot() #using PyPlot backend

# Set a random seed for reproducible behaviour
rng = StableRNG(1111)

#Generate the training data

function lotka!(du, u, p, t)
    α, β, γ, δ = p
    du[1] = α * u[1] - β * u[2] * u[1]
    du[2] = γ * u[1] * u[2] - δ * u[2]
end

# Define the experimental parameter
tspan = (0.0, 5.0)
u0 = 5.0f0 * rand(rng, 2)
p_ = [1.3, 0.9, 0.8, 1.8]           #truth parameters
prob = ODEProblem(lotka!, u0, tspan, p_)
solution = solve(prob, Vern7(), abstol = 1e-12, reltol = 1e-12, saveat = 0.25)

# Add noise in terms of the mean
X = Array(solution)
t = solution.t

x̄ = mean(X, dims = 2)
noise_magnitude = 5e-3
Xₙ = X .+ (noise_magnitude * x̄) .* randn(rng, eltype(X), size(X))

my_plot = plot(solution, alpha = 0.75, color = :black, label = ["True Data" nothing])
scatter!(t, transpose(Xₙ), color = :red, label = ["Noisy Data" nothing])
display(my_plot)

# Define UDE 
rbf(x) = exp.(-(x .^ 2))

# Multilayer FeedForward
const U = Lux.Chain(Lux.Dense(2, 5, rbf), Lux.Dense(5, 5, rbf), Lux.Dense(5, 5, rbf),
    Lux.Dense(5, 2))
# Get the initial parameters and state variables of the model
p, st = Lux.setup(rng, U)
const _st = st

# Define the hybrid model
function ude_dynamics!(du, u, p, t, p_true)
    û = U(u, p, _st)[1] # Network prediction
    du[1] = p_true[1] * u[1] + û[1]
    du[2] = -p_true[4] * u[2] + û[2]
end

# Closure with the known parameter
nn_dynamics!(du, u, p, t) = ude_dynamics!(du, u, p, t, p_)
# Define the problem
prob_nn = ODEProblem(nn_dynamics!, Xₙ[:, 1], tspan, p)


# setup the training loop
function predict(θ, X = Xₙ[:, 1], T = t)
    _prob = remake(prob_nn, u0 = X, tspan = (T[1], T[end]), p = θ)
    Array(solve(_prob, Vern7(), saveat = T,
        abstol = 1e-6, reltol = 1e-6,
        sensealg = QuadratureAdjoint(autojacvec = ReverseDiffVJP(true))))
end

# Define the loss function
function loss(θ)
    X̂ = predict(θ)
    mean(abs2, Xₙ .- X̂)
end

# define callback
losses = Float64[]

callback = function (state, l)
    push!(losses, l)
    if length(losses) % 50 == 0
        println("Current loss after $(length(losses)) iterations: $(losses[end])")
    end
    return false
end

# Define the optimization problem
adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)

optprob = Optimization.OptimizationProblem(optf, ComponentVector{Float64}(p))

# Solve the optimization problem first with the Adam optimizer
res1 = Optimization.solve(
    optprob, OptimizationOptimisers.Adam(), callback = callback, maxiters = 5000)
println("Training loss after $(length(losses)) iterations: $(losses[end])\n")

# Solve the optimization problem with the LBFGS optimizer
optprob2 = Optimization.OptimizationProblem(optf, res1.u)
res2 = Optimization.solve(
    optprob2, LBFGS(linesearch = BackTracking()), callback = callback, maxiters = 1000)
println("Final training loss after $(length(losses)) iterations: $(losses[end])")

# Rename the best candidate
p_trained = res2.u

@show typeof(p_trained)
@show p_trained

#save the model
jldsave("model.jld", ps=p_trained, st=st)

#load the model
p_trained_loaded, st_loaded = jldopen("model.jld", "r") do file
    file["ps"], file["st"]
end

println(" ")
@show typeof(p_trained_loaded)
@show p_trained_loaded

println(" ")
#compare the loaded model with the trained model
@show p_trained == p_trained_loaded

# Plot the results
# Plot the losses
#pl_losses = plot(1:5000, losses[1:5000], yaxis = :log10, xaxis = :log10,
#    xlabel = "Iterations", ylabel = "Loss", label = "ADAM", color = :blue)
#plot!(5001:length(losses), losses[5001:end], yaxis = :log10, xaxis = :log10,
#    xlabel = "Iterations", ylabel = "Loss", label = "LBFGS", color = :red)

#display(pl_losses)

# Plot the true and predicted data
## Analysis of the trained network
# Plot the data and the approximation
ts = first(solution.t):(mean(diff(solution.t)) / 2):last(solution.t)
X̂ = predict(p_trained, Xₙ[:, 1], ts)
# Trained on noisy data vs real solution
pl_trajectory = plot(ts, transpose(X̂), xlabel = "t", ylabel = "x(t), y(t)", color = :red,
    label = ["UDE Approximation" nothing])
scatter!(solution.t, transpose(Xₙ), color = :black, label = ["Measurements" nothing])

display(pl_trajectory)

# plot unknown terms
# Ideal unknown interactions of the predictor
Ȳ = [-p_[2] * (X̂[1, :] .* X̂[2, :])'; p_[3] * (X̂[1, :] .* X̂[2, :])']
# Neural network guess
Ŷ = U(X̂, p_trained, st)[1]

pl_reconstruction = plot(ts, transpose(Ŷ), xlabel = "t", ylabel = "U(x,y)", color = :red,
    label = ["UDE Approximation" nothing])
plot!(ts, transpose(Ȳ), color = :black, label = ["True Interaction" nothing])

display(pl_reconstruction)

gui()

readline()

println("Done!")