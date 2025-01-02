using DifferentialEquations, Debugger

function my_ode!(du, u, p, t)
    du[1] = -p[1] * u[1]
end

u0 = [1.0]            # Initial condition
p = [0.5]             # Parameters
tspan = (0.0, 10.0)   # Time span

# Define the ODE problem
prob = ODEProblem(my_ode!, u0, tspan, p)

# Debug the solver
@enter sol = solve(prob, Tsit5())

#after debugger
println("sol: ", sol)
