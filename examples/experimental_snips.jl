using Revise
using DifferentialEquations
using SciMLSensitivity
using Zygote
using ComponentArrays
using Optimization, OptimizationOptimisers

# Define the shallow water equations (non-mutating version)
function shallow_water_eq(u, p, t)
    # Extract parameters from ComponentArray
    zb, n, Q = p.zb[1], p.n, p.Q
    h, hu, hv = u
    g = 9.81

    Zygote.ignore() do
        println("zb = ", zb)
        println("n = ", n)
        println("Q = ", Q)
        println("sum(Q) = ", sum(Q))
        println("h = ", h)
        println("hu = ", hu)
        println("hv = ", hv)
    end

    # Return derivatives without mutation
    [
        sum(Q) - hu - hv - zb,  # Example using Q as a spatial inflow array
        -g * h - n * hu,
        -g * h - n * hv
    ]
end

# Define the loss function (non-mutating version)
function compute_loss(p, u0, tspan, observed_data)
    # Create ODEProblem with non-mutating function
    prob = ODEProblem(shallow_water_eq, u0, tspan, p)
    
    # Solve the ODE
    sol = solve(prob, Tsit5(), saveat=0.1)
    
    # Loss: Mean squared error between simulation and observation
    simulated = sol[end]
    sum((simulated .- observed_data).^2)
end

# Setup optimization (non-mutating version)
function optimize_parameters(u0, tspan, observed_data, p_init, active_params; maxiters=100)
    # Loss function for optimization
    function opt_loss(θ,p)
        # Create new ComponentArray based on active parameters
        if active_params == [:zb]
            p_new = ComponentArray(zb=θ, n=p_init.n, Q=p_init.Q)
        elseif active_params == [:n]
            p_new = ComponentArray(zb=p_init.zb, n=θ, Q=p_init.Q)
        elseif active_params == [:Q]
            p_new = ComponentArray(zb=p_init.zb, n=p_init.n, Q=θ)
        elseif active_params == [:zb, :Q]
            p_new = ComponentArray(zb=θ[1:1], n=p_init.n, Q=θ[2:end])
        else
            error("Unsupported parameter combination")
        end
        
        compute_loss(p_new, u0, tspan, observed_data)
    end

    # Initial parameter values (flattened)
    θ0 = vcat([param == :Q ? p_init.Q : [p_init[param]] for param in active_params]...)

    # Create optimization problem
    optf = OptimizationFunction(opt_loss, Optimization.AutoZygote())
    optprob = OptimizationProblem(optf, θ0)

    # Solve optimization problem
    solve(optprob, Adam(0.01), maxiters=maxiters)
end

# Test the optimization
u0 = [1.0, 0.1, 0.1]  # Initial conditions: [h, hu, hv]
tspan = (0.0, 10.0)   # Time span
observed_data = [0.9, 0.05, 0.05]  # Example observed data

# Create initial parameters
p_init = ComponentArray(zb=0.5, n=0.03, Q=[0.5, 0.6, 0.7])
active_params = [:zb, :Q]  # Optimize both zb and Q

# Optimize parameters
p_optimized = optimize_parameters(u0, tspan, observed_data, p_init, active_params)
println("Optimized Parameters: ", p_optimized)