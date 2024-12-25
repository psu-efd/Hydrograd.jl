using Revise
using DifferentialEquations
using SciMLSensitivity
using Zygote
using ComponentArrays  # Add this
using Optimization, OptimizationOptimisers  # Add these

# Define the shallow water equations
function shallow_water_eq!(du, u, p, t)
    # Extract parameters from ComponentArray
    zb, n, Q = p.zb, p.n, p.Q
    h, hu, hv = u
    g = 9.81

    # Compute derivatives
    du[1] = Q - hu - hv - zb
    du[2] = -g * h - n * hu
    du[3] = -g * h - n * hv
end

# Define the loss function
function compute_loss(p, u0, tspan, observed_data)
    # Create ODEProblem
    prob = ODEProblem(shallow_water_eq!, u0, tspan, p)
    
    # Solve the ODE
    sol = solve(prob, Tsit5(), saveat=0.1)
    
    # Loss: Mean squared error between simulation and observation
    simulated = sol[end]
    return sum((simulated .- observed_data).^2)
end

# Example Usage
u0 = [1.0, 0.1, 0.1]  # Initial conditions: [h, hu, hv]
tspan = (0.0, 10.0)    # Time span
observed_data = [0.9, 0.05, 0.05]  # Example observed data

# Create ComponentArray for parameters
p = ComponentArray(zb=0.5, n=0.03, Q=1.0)
active_params = [:zb, :Q]  # Parameters to optimize

# Compute gradients
function compute_gradient(u0, tspan, observed_data, p, active_params)
    # Create function that only depends on active parameters
    function loss_active(active_vals)
        # Create new parameter set without mutation
        new_params = Dict(zip(active_params, active_vals))
        p_new = ComponentArray(
            zb = haskey(new_params, :zb) ? new_params[:zb] : p.zb,
            n = haskey(new_params, :n) ? new_params[:n] : p.n,
            Q = haskey(new_params, :Q) ? new_params[:Q] : p.Q
        )
        return compute_loss(p_new, u0, tspan, observed_data)
    end
    
    # Extract active parameter values
    active_vals = [p[param] for param in active_params]
    
    # Compute gradient
    loss, back = Zygote.pullback(loss_active, active_vals)
    gradient = back(1.0)[1]
    
    return loss, gradient
end

# Test the gradient computation
loss, gradient = compute_gradient(u0, tspan, observed_data, p, active_params)
println("Loss: ", loss)
println("Gradient: ", gradient)

# Setup optimization
function optimize_parameters(u0, tspan, observed_data, p_init, active_params; maxiters=100)
    # Loss function for optimization
    function opt_loss(θ, p)  # Add p argument even if unused
        # Create parameter set from optimization variables
        new_params = Dict(zip(active_params, θ))
        p_new = ComponentArray(
            zb = haskey(new_params, :zb) ? new_params[:zb] : p_init.zb,
            n = haskey(new_params, :n) ? new_params[:n] : p_init.n,
            Q = haskey(new_params, :Q) ? new_params[:Q] : p_init.Q
        )
        return compute_loss(p_new, u0, tspan, observed_data)
    end

    # Initial values for optimization
    θ0 = [p_init[param] for param in active_params]

    # Create optimization problem
    optf = OptimizationFunction(opt_loss, Optimization.AutoZygote())
    optprob = OptimizationProblem(optf, θ0)  # No parameters needed

    # Solve optimization problem
    sol = solve(optprob, Adam(0.01), maxiters=maxiters)
    
    return sol
end

# Test the optimization
p_init = ComponentArray(zb=0.0, n=0.03, Q=0.5)  # Initial guess
sol = optimize_parameters(u0, tspan, observed_data, p_init, active_params)
println("Optimized parameters: ", sol.u)