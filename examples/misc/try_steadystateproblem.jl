using DifferentialEquations
using OrdinaryDiffEq
using SteadyStateDiffEq
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using ReverseDiff

SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true

# Define the system of equations
function f(du, u, p, t)
    # Bound the values to prevent overflow
    u1_bounded = clamp(u[1], -100.0, 100.0)
    u2_bounded = clamp(u[2], -100.0, 100.0)
    
    du[1] = -p[1]*u1_bounded + u2_bounded^2 - 1.0
    du[2] = -p[2]*u2_bounded + sin(u1_bounded) - 2.0  # sin now has bounded input
end

# Initial conditions and parameters
u0 = [0.5, 0.5]
p_true = [1.0, 1.0]
p_init = [0.4, 0.4]

# Generate synthetic data
steady_prob = SteadyStateProblem(f, u0, p_true)
sol_true = solve(steady_prob, DynamicSS(Tsit5()))
u_true = sol_true.u

# Define loss function with error handling
function loss_function(p, _)
    steady_prob = SteadyStateProblem(f, u0, p)
    
    # Add solver options for better stability
    sol = try
        solve(steady_prob, DynamicSS(Tsit5()), 
              sensealg=SteadyStateAdjoint(autojacvec=true),
              abstol=1e-8, reltol=1e-6,
              maxiters=1000,
              verbose=false)
    catch e
        @warn "Solver failed" p
        return convert(eltype(p), Inf)  # Return large loss for failed solutions
    end

    @show sol
    
    # Check if solution is valid
    #if sol.retcode == SciMLBase.ReturnCode.Success
    #    return convert(eltype(p), Inf)
    #end
    
    u_pred = sol.u
    return sum(abs2, u_pred .- u_true)
end

# Setup optimization problem with bounds
adtype = Optimization.AutoZygote()
optf = OptimizationFunction(loss_function, adtype)
optprob = OptimizationProblem(optf, p_init, 
                             lb=[0.0, 0.0],
                             ub=[10.0, 10.0])

# Solve with more robust settings
sol = solve(optprob, Optim.Fminbox(BFGS()), 
           maxiters=100,
           allow_f_increases=true,
           f_abstol=1e-8)

# Print results
println("True parameters: ", p_true)
println("Found parameters: ", sol.u)
println("Final loss: ", sol.minimum)

# Verify solution
steady_prob_final = SteadyStateProblem(f, u0, sol.u)
sol_final = solve(steady_prob_final, DynamicSS(Tsit5()))
println("True steady state: ", u_true)
println("Final steady state: ", sol_final.u)
