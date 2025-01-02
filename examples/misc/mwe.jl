using SciMLSensitivity

using DifferentialEquations

using Optimization
using OptimizationOptimisers

using Plots

#AD engines
using Zygote
using ForwardDiff
using ReverseDiff

      #directory to save to (the same directory as the main file)
save_path = dirname(@__FILE__)

function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p

    #test whether if-else breaks AD
    if α < 1.0
        α = 1.2
    else
        α = 1.5
    end

    du[1] = dx = α * x - β * x * y
    du[2] = dy = -δ * y + γ * x * y

    #println("du/dp = ", ForwardDiff.partials(du))
end

# Initial condition
u0 = [1.0, 1.0]

# Simulation interval and intermediary points
tspan = (0.0, 4.0)
tsteps = 0.0:0.1:4.0

# LV equation parameter. p = [α, β, δ, γ]
p = [1.5, 1.0, 3.0, 1.0]

# Setup the ODE problem, then solve
prob = ODEProblem(lotka_volterra!, u0, tspan, p)
sol_true = solve(prob, Tsit5(), p = p, saveat = tsteps)

# Plot the solution
#println("plotting true solution...")
#plot(sol_true)
#savefig(joinpath(save_path, "LV_ode.png"))


function loss(p)
    sol = solve(prob, Tsit5(), p = p, saveat = tsteps)
    #sol = solve(prob, TRBDF2(), p = p, saveat = tsteps)
    pred = Array(sol)
    loss = sum(abs2, pred .- Array(sol_true))
    #println(pred[end,:])
    #println(Array(sol_true)[end,:])

    Zygote.ignore() do
        #pred_real = ForwardDiff.value.(pred)
        #println("size of pred_real", size(pred_real))
        #println("size of sol", size(sol.u))
        #println("sol.u[1] ", ForwardDiff.value.(sol.u[1]))
        #println("sol.u[end] ", ForwardDiff.partials.(sol.u[end]))
        #println("sol.t ", sol.t)
        #println("size of sol_true", size(Array(sol_true)))
        #println("loss = ", loss)

        # display(loss)
        # plt = plot(pred_real', 
        #     label="Predicted",
        #     #ylim=(0, 10),
        #     xlabel="t",
        #     ylabel="u(t)",
        #     linestyle=:solid)

        # plot!(Array(sol_true)', label="True", linestyle=:dash)

        # display(plt)
    end

    return loss
end

#set inital guess for inversion
#p_init = [1.4, 1.0, 2.0, 1.0]
p_init = [0.0, 0.0, 0.0, 0.0]

# Add this before gradient computation
#@code_warntype(loss(p_init))

SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true
#grad = ForwardDiff.gradient(loss, p_init)
#grad = Zygote.gradient(loss, p_init)
grad = ReverseDiff.gradient(loss, p_init)
#jac = Zygote.jacobian(predict, ps)
#jac = ReverseDiff.jacobian(predict, ps)
@show grad
#plot(jac)
#readline()
throw("stop here")

callback = function (state, l)
    display(l)
    pred = solve(prob, Tsit5(), p = state.u, saveat = tsteps)
    #plt = plot(pred, ylim = (0, 6))
    plt = plot(pred, 
         label="Predicted",
         ylim=(0, 10),
         xlabel="t",
         ylabel="u(t)",
         linestyle=:solid)

    plot!(sol_true, label="True", linestyle=:dash)

    display(plt)

    println("loss = ", l)
    println("p = ", state.u)

    println("Press Enter to continue...")
    readline()

    # Tell Optimization.solve to not halt the optimization. If return true, then
    # optimization stops.
    return false
end



adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, p_init)

result_ode = Optimization.solve(optprob, PolyOpt(), callback = callback,  maxiters = 500)
#result_ode = Optimization.solve(optprob, Adam(0.01), callback = callback,  maxiters = 500)