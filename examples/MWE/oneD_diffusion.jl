using DifferentialEquations
using Optimization
using OptimizationOptimisers
using Plots
using Zygote

# Define the problem parameters
nx = 50                  # Number of spatial points
L = 1.0                 # Domain length
dx = L/nx
x = range(0, L, nx)
tspan = (0.0, 1.0)

# True diffusion coefficient (we'll try to recover this)
D_true = 0.1

# Create synthetic data
function diffusion_equation!(du, u, p, t)
    D = p[1]  # Diffusion coefficient
    # Interior points
    for i in 2:nx-1
        du[i] = D * (u[i+1] - 2u[i] + u[i-1]) / dx^2
    end
    # Boundary conditions (fixed at 0)
    du[1] = 0.0
    du[nx] = 0.0
end

# Initial condition (Gaussian pulse)
u0 = exp.(-100 * (x .- 0.5).^2)

# Create and solve the problem with true diffusion coefficient
prob_true = ODEProblem(diffusion_equation!, u0, tspan, [D_true])
sol_true = solve(prob_true, Tsit5(), saveat=0.1)

# Generate synthetic data with some noise
data = Array(sol_true) + 0.01 * randn(size(Array(sol_true)))

# Define the loss function using only the final time step
function loss(p, _)
    prob = ODEProblem(diffusion_equation!, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=sol_true.t)
    pred = Array(sol)
    return sum(abs2, pred[end,:] .- Array(sol_true)[end,:])  # Compare only final time step
end

# Initial guess for diffusion coefficient
D_init = [0.2]  # Starting with wrong guess

# Define the optimization problem with explicit in-place specification
optf = OptimizationFunction(loss, Optimization.AutoZygote(); grad=false)
optprob = OptimizationProblem(optf, D_init)

# Callback to track progress
losses = []
Ds = []
callback = function (p, l)
    push!(losses, l)
    push!(Ds, p[1])
    println("D = $(p[1]), Loss = $l")
    return false
end

# Solve the optimization problem
result = Optimization.solve(optprob, Adam(0.01), callback=callback, maxiters=100)

# Plot results
p1 = plot(x, Array(sol_true)[end, :], label="True solution", title="Final State")
plot!(p1, x, Array(solve(ODEProblem(diffusion_equation!, u0, tspan, [result.u[1]]), Tsit5(), saveat=0.1))[end, :], 
      label="Recovered solution")

p2 = plot(losses, yscale=:log10, xlabel="Iteration", ylabel="Loss", label="")

p3 = plot(Ds, ylabel="Diffusion coefficient", xlabel="Iteration", label="",
          title="D_true = $D_true, D_recovered = $(round(result.u[1], digits=4))")

plot(p1, p2, p3, layout=(3,1), size=(800,1000))