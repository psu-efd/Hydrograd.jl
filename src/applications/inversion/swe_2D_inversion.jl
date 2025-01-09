
#perform inversion for 2D SWE
using Dates

#using Debugger, JuliaInterpreter

#main function for the inversion
function swe_2D_inversion(ode_f, Q0, params_vector, swe_extra_params)

    #unpack the extra parameters
    active_param_name = swe_extra_params.active_param_name
    settings = swe_extra_params.settings
    my_mesh_2D = swe_extra_params.my_mesh_2D
    swe_2D_constants = swe_extra_params.swe_2D_constants
    case_path = swe_extra_params.case_path

    nodeCoordinates = swe_extra_params.nodeCoordinates

    if settings.bVerbose
        println("       inversion parameter name = ", active_param_name)
    end

    #some sanity check
    if settings.inversion_settings.bInversion_slope_loss && active_param_name !== "zb"
        throw(ArgumentError("The slope regularization is only applicable to the bed elevation parameter (zb). No inversion is performed. Please check the inversion settings."))
    end

    #open the forward simulation result (as the ground truth)
    sol_truth = JSON3.read(open(joinpath(case_path, settings.inversion_settings.inversion_truth_file_name)), Dict{String,Vector{Float64}})

    if settings.bVerbose
        #@show typeof(sol_truth)
        #@show length(sol_truth)
        #@show sol_truth
    end

    WSE_truth = vec(sol_truth["wse_truth"])
    h_truth = vec(sol_truth["h_truth"])
    u_truth = vec(sol_truth["u_truth"])
    v_truth = vec(sol_truth["v_truth"])
    zb_cell_truth = vec(sol_truth["zb_cell_truth"])
    ManningN_cell_truth = vec(sol_truth["ManningN_cells_truth"])
    inlet_discharges_truth = vec(sol_truth["inlet_discharges_truth"])

    println("   Loading inversion data ...")

    #combine the truth data into a dictionary
    observed_data = Dict{String,Vector{Float64}}("WSE_truth" => WSE_truth, "h_truth" => h_truth, "u_truth" => u_truth, "v_truth" => v_truth,
        "zb_cell_truth" => zb_cell_truth, "ManningN_cell_truth" => ManningN_cell_truth, "inlet_discharges_truth" => inlet_discharges_truth)

    #debug AD correctness
    #debug start
    #debug_AD(ode_f, Q0, swe_2D_constants, params_vector, swe_extra_params)
    #return
    #debug end

    Zygote.ignore() do
        if settings.bVerbose
            @show typeof(Q0)
            @show typeof(swe_2D_constants.tspan)
            @show typeof(observed_data)
            @show typeof(params_vector)

            @show size(params_vector)
            @show params_vector
        end
    end

    #perform the inversion
    println("   Performing inversion ...\n")
    sol, ITER, LOSS, PRED, PARS = optimize_parameters_inversion(ode_f, Q0, swe_2D_constants.tspan, params_vector, settings, my_mesh_2D, swe_2D_constants, observed_data, active_param_name, case_path)

    #save the inversion results
    jldsave(joinpath(case_path, settings.inversion_settings.save_file_name); ITER, LOSS, PRED, PARS)

    #process the inversion results
    println("   Post-processing inversion results ...")

    #process inversion results
    Hydrograd.postprocess_inversion_results_swe_2D(settings, my_mesh_2D, nodeCoordinates, zb_cell_truth, h_truth, u_truth, v_truth, WSE_truth, case_path)

end


# Define the loss function
function compute_loss_inversion(ode_f::ODEFunction, Q0::Matrix{T1}, tspan::Tuple{T1,T1}, p::AbstractVector{T2}, settings::ControlSettings, 
    my_mesh_2D::mesh_2D, swe_2D_constants::swe_2D_consts, observed_data::Dict{String,Vector{T1}}, active_param_name::String, data_type::DataType) where {T1<:Real,T2<:Real}

    # Solve the ODE (forward pass)

    # Create ODEProblem        
    prob = ODEProblem(ode_f, Q0, tspan, p)

    # Create ODE solver
    if settings.inversion_settings.ode_solver == "Tsit5()"
        ode_solver = Tsit5()
    else
        error("Not implemented yet")
    end

    # Create sensealg for the ODE solver
    ode_solver_sensealg = nothing  #default sensealg is nothing
    #For sensealg argument in ODE solve function: See https://docs.sciml.ai/SciMLSensitivity/dev/manual/differential_equation_sensitivities/ for the choice of AD type (sensealg)
    #If not specified, the default is a smart polyalgorithm used to automatically determine the most appropriate method for a given equation. It is defined in 
    #SciMLSensitivity.jl/src/concrete_solve.jl 
    #default_sensealg = if p !== SciMLBase.NullParameters() &&
    #                  !(eltype(u0) <: ForwardDiff.Dual) &&
    #                  !(eltype(p) <: ForwardDiff.Dual) &&
    #                  !(eltype(u0) <: Complex) &&
    #                  !(eltype(p) <: Complex) &&
    #                  length(u0) + length(tunables) <= 100
    #    ForwardDiffSensitivity()    #for small problems, it uses ForwardDiffSensitivity()
    #    ...
    #Some options are:
    # 1. BacksolveAdjoint(autojacvec=ZygoteVJP())       #working
    # 2. InterpolatingAdjoint(autojacvec=ZygoteVJP())   #working
    # 3. ReverseDiffAdjoint(autojacvec=ZygoteVJP())     #not working
    # 4. TrackerAdjoint(autojacvec=ZygoteVJP())        #not tested
    # 5. ZygoteVJP()                                  #working
    # 6. TrackerVJP()                                 #not tested
    # 7. ReverseDiffVJP()                             #not tested
    # 8. InterpolatingVJP()                          #not tested
    # 9. BacksolveVJP()                              #not tested
    # 10. Default (unspecified): working, but no control on sensealg
    if settings.inversion_settings.ode_solver_sensealg == "ZygoteVJP()"
        ode_solver_sensealg = ZygoteVJP()
    elseif settings.inversion_settings.ode_solver_sensealg == "InterpolatingAdjoint(autojacvec=ZygoteVJP())"
        ode_solver_sensealg = InterpolatingAdjoint(autojacvec=ZygoteVJP())
    elseif settings.inversion_settings.ode_solver_sensealg == "AutoForwardDiff()"
        ode_solver_sensealg = ForwardDiffSensitivity()
    elseif settings.inversion_settings.ode_solver_sensealg == "BacksolveAdjoint(autojacvec=ZygoteVJP())"
        ode_solver_sensealg = BacksolveAdjoint(autojacvec=ZygoteVJP())
    else
        error("Not implemented yet")
    end

    #define the time for saving the results for the ODE solver
    dt_save = (swe_2D_constants.tspan[2] - swe_2D_constants.tspan[1]) / settings.inversion_settings.ode_solver_nSave
    t_save = swe_2D_constants.tspan[1]:dt_save:swe_2D_constants.tspan[2]

    # Solve the ODE
    pred = solve(prob, ode_solver, adaptive=settings.inversion_settings.ode_solver_adaptive, dt=swe_2D_constants.dt, saveat=t_save, sensealg=ode_solver_sensealg)

    Zygote.ignore() do
        if settings.bVerbose
            #@show pred.retcode
            #@show typeof(pred)
            #@show size(pred)
            #@show pred
        end
    end

    #compute the loss
    # Ensure type stability in loss computation
    loss_total = zero(data_type)
    loss_pred = zero(data_type)
    loss_pred_WSE = zero(data_type)
    loss_pred_uv = zero(data_type)
    loss_bound = zero(data_type)
    loss_slope = zero(data_type)

    if pred.retcode == SciMLBase.ReturnCode.Success
        WSE_truth = Vector{Float64}(observed_data["WSE_truth"])
        h_truth = Vector{Float64}(observed_data["h_truth"])
        u_truth = Vector{Float64}(observed_data["u_truth"])
        v_truth = Vector{Float64}(observed_data["v_truth"])

        #Get the bed elevation at cells: if active_param_name == "zb", then use the initial values; otherwise, use the truth data
        if active_param_name == "zb"
            zb_cells_temp = p
        else
            zb_cells_temp = observed_data["zb_cell_truth"]
        end

        l = pred[:, 1, end] .+ zb_cells_temp .- WSE_truth  #loss = free surface elevation mismatch

        #loss for free surface elevation mismatch
        if settings.inversion_settings.bInversion_WSE_loss
            loss_pred_WSE = sum(abs2, l)
        end

        #loss for velocity mismatch
        # Add small epsilon to prevent division by zero
        ϵ = sqrt(eps(data_type))

        if settings.inversion_settings.bInversion_uv_loss      #if also include u in the loss 
            l_u = pred[:, 2, end] ./ (pred[:, 1, end] .+ ϵ) .- u_truth
            l_v = pred[:, 3, end] ./ (pred[:, 1, end] .+ ϵ) .- v_truth

            loss_pred_uv = sum(abs2, l_u) + sum(abs2, l_v)
        end

        #combined loss due to free surface elevation mismatch and velocity mismatch
        loss_pred = loss_pred_WSE + loss_pred_uv

        #loss for parameter bound regularization
        if settings.inversion_settings.bInversion_bound_loss
            loss_bound = compute_bound_loss(p, settings.inversion_settings.parameter_value_lower_bound, settings.inversion_settings.parameter_value_upper_bound)
        end

        #loss for bed slope regularization
        if settings.inversion_settings.bInversion_slope_loss    #if bed slope is included in the loss 
            loss_slope = calc_slope_loss(zb_cells_temp, settings, my_mesh_2D)
        end

        #combined loss due to free surface elevation mismatch, velocity mismatch, and bed slope regularization
        loss_total = loss_pred + loss_bound + loss_slope
    else
        loss_total = convert(data_type, Inf)
    end

    Zygote.ignore() do
        if settings.bVerbose
            #@show loss_total
            #@show loss_pred
            #@show loss_pred_WSE
            #@show loss_pred_uv
            #@show loss_slope
        end
    end

    return loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_bound, loss_slope, pred
end

function optimize_parameters_inversion(ode_f, Q0, tspan, p_init, settings, my_mesh_2D, swe_2D_constants, observed_data, active_param_name, case_path)
    #start the timer for the inversion
    inversion_start_time = now()  # Current date and time

    # Loss function for optimization
    function opt_loss(θ, p)  # Add p argument even if unused

        data_type = eltype(θ)

        Zygote.ignore() do
            #println("data_type = ", data_type)
            #println("active_range = ", active_range)
            #println("param_ranges = ", param_ranges)
            #println("θ = ", θ)
            #println("p_init = ", p_init)
        end

        loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_bound, loss_slope, pred = compute_loss_inversion(ode_f, Q0, tspan, θ, settings,
            my_mesh_2D, swe_2D_constants, observed_data, active_param_name, data_type)

        return loss_total
    end

    # Define AD type choice for optimization's gradient computation
    #The following is from SciMLSensitivity documentation regarding the choice of AD 
    #  AutoForwardDiff(): The fastest choice for small optimizations
    #  AutoReverseDiff(compile=false): A fast choice for large scalar optimizations
    #  AutoTracker(): Like ReverseDiff but GPU-compatible
    #  AutoZygote(): The fastest choice for non-mutating array-based (BLAS) functions
    #  AutoFiniteDiff(): Finite differencing, not optimal but always applicable
    #  AutoModelingToolkit(): The fastest choice for large scalar optimizations
    #  AutoEnzyme(): Highly performant AD choice for type stable and optimized code

    adtype = nothing
    if settings.inversion_settings.inversion_sensealg == "AutoZygote()"
        adtype = Optimization.AutoZygote()
    elseif settings.inversion_settings.inversion_sensealg == "AutoReverseDiff()"
        adtype = Optimization.AutoReverseDiff(compile=false)
    elseif settings.inversion_settings.inversion_sensealg == "AutoForwardDiff()"
        adtype = Optimization.AutoForwardDiff()
    elseif settings.inversion_settings.inversion_sensealg == "AutoReverseDiff()"
        adtype = Optimization.AutoReverseDiff(compile=false)
    else
        println("       settings.inversion_settings.inversion_sensealg = ", settings.inversion_settings.inversion_sensealg)
        throw(ArgumentError("Invalid sensealg choice. Supported sensealg: AutoZygote(), AutoReverseDiff(), AutoForwardDiff(), AutoFiniteDiff(), AutoModelingToolkit(), AutoEnzyme(). No inversion is performed."))
    end

    # Define the optimization function and problem
    #From SciMLSensitivity documentation: https://docs.sciml.ai/Optimization/stable/API/optimization_function/
    # OptimizationFunction{iip}(f, adtype::AbstractADType = NoAD();
    #                       grad = nothing, hess = nothing, hv = nothing,
    #                       cons = nothing, cons_j = nothing, cons_jvp = nothing,
    #                       cons_vjp = nothing, cons_h = nothing,
    #                       hess_prototype = nothing,
    #                       cons_jac_prototype = nothing,
    #                       cons_hess_prototype = nothing,
    #                       observed = __has_observed(f) ? f.observed : DEFAULT_OBSERVED_NO_TIME,
    #                       lag_h = nothing,
    #                       hess_colorvec = __has_colorvec(f) ? f.colorvec : nothing,
    #                       cons_jac_colorvec = __has_colorvec(f) ? f.colorvec : nothing,
    #                       cons_hess_colorvec = __has_colorvec(f) ? f.colorvec : nothing,
    #                       lag_hess_colorvec = nothing,
    #                       sys = __has_sys(f) ? f.sys : nothing)
    optf = OptimizationFunction(opt_loss, adtype)

    #From SciMLSensitivity documentation: https://docs.sciml.ai/Optimization/stable/API/optimization_problem/
    #OptimizationProblem{iip}(f, u0, p = SciMLBase.NullParameters(),;
    #                          lb = nothing,
    #                          ub = nothing,
    #                          lcons = nothing,
    #                          ucons = nothing,
    #                          sense = nothing,
    #                          kwargs...)
    optprob = OptimizationProblem(optf, p_init)  # No parameters needed

    # Define the bounds for the parameter (only applicable for some optimizers which support lb and ub)
    # Note: lb and ub are not used for the optimizer currently
    #lb_p = zeros(my_mesh_2D.numOfCells)
    #lb_p .= -0.1
    #ub_p = zeros(my_mesh_2D.numOfCells)
    #ub_p .= 0.3


    #define the accumulators for the inversion results
    ITER = []        #iteration number accumulator
    LOSS = []        # Loss accumulator
    PRED = []        # prediction accumulator
    PARS = []        # parameters accumulator

    # Define the callback to handle both vector and OptimizationState inputs
    #callback = function (state_or_θ, loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope, pred)
    #    θ = state_or_θ isa Optimization.OptimizationState ? state_or_θ.u : state_or_θ
    #    println("Loss: ", loss_total)
    #    println("Parameters: ", θ)
    #    return false  # continue optimization
    #end

    #callback function to observe training
    callback = function (optimizer_state, loss_total)

        inversion_current_time = now()  # Current date and time
        inversion_elapsed_time = inversion_current_time - inversion_start_time
        inversion_elapsed_seconds = Millisecond(inversion_elapsed_time).value / 1000

        #reset the inversion start time
        inversion_start_time = inversion_current_time

        θ = optimizer_state.u

        iter_number = length(LOSS) + 1  #get the inversion iteration number 

        Zygote.ignore() do
            println("       iter_number = ", iter_number, ", loss_total = ", loss_total, ", inversion_elapsed_seconds = ", inversion_elapsed_seconds)

            #save the inversion results
            if iter_number % settings.inversion_settings.save_frequency == 0

                data_type = eltype(θ)

                loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_bound, loss_slope, prediction = compute_loss_inversion(ode_f, Q0, tspan, θ, settings,
                        my_mesh_2D, swe_2D_constants, observed_data, active_param_name, data_type)

                append!(ITER, iter_number)

                append!(PRED, [prediction[:, :, end]])

                append!(LOSS, [[loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_bound, loss_slope]])

                append!(PARS, [θ])
            end

            #checkpoint the inversion results (in case the inversion is interrupted)
            if settings.inversion_settings.save_checkpoint && iter_number % settings.inversion_settings.checkpoint_frequency == 0
                jldsave(joinpath(case_path, "checkpoint_inversion_iter_$(iter_number).jld2"); ITER, LOSS, PRED, PARS)
            end

        end



        #if l > 1e-9
        #    false
        #else 
        #    true   #force the optimizer to stop 
        #end

        false
    end

    sol = nothing

    #loop over the optimizers
    for (iOptimizer, optimizer_choice) in enumerate(settings.inversion_settings.inversion_optimizers)
        println("   iOptimizer: ", iOptimizer, " out of ", length(settings.inversion_settings.inversion_optimizers))
        println("       optimizer_choice: ", optimizer_choice)
        println("       iterations: ", settings.inversion_settings.inversion_max_iterations[iOptimizer])
        println("       learning rate: ", settings.inversion_settings.inversion_learning_rates[iOptimizer])
        println("       abs_tol: ", settings.inversion_settings.inversion_abs_tols[iOptimizer])
        println("       rel_tol: ", settings.inversion_settings.inversion_rel_tols[iOptimizer])

        max_iterations = settings.inversion_settings.inversion_max_iterations[iOptimizer]

        if max_iterations <= 0
            println("       max_iterations <= 0. No inversion is performed with this optimizer.")
            continue
        end

        #create the optimizer
        if optimizer_choice == "Adam"
            optimizer = OptimizationOptimisers.Adam(settings.inversion_settings.inversion_learning_rates[iOptimizer])
        elseif optimizer_choice == "LBFGS"
            optimizer = OptimizationOptimJL.LBFGS(linesearch=LineSearches.BackTracking())
        else
            throw(ArgumentError("Invalid optimizer choice. Supported optimizers: Adam, LBFGS. No inversion is performed."))
        end

        # Solve optimization problem
        # From SciMLSensitivity documentation: https://docs.sciml.ai/Optimization/stable/API/optimization_solution/
        # Returned optimization solution Fields:
        #   u: the representation of the optimization's solution.
        #   cache::AbstractOptimizationCache: the optimization cache` that was solved.
        #   alg: the algorithm type used by the solver.
        #   objective: Objective value of the solution
        #   retcode: the return code from the solver. Used to determine whether the solver solved successfully or whether 
        #            it exited due to an error. For more details, see the return code documentation.
        #   original: if the solver is wrapped from a external solver, e.g. Optim.jl, then this is the original return from said solver library.
        #   stats: statistics of the solver, such as the number of function evaluations required.

        if optimizer_choice == "Adam" 
            sol = solve(optprob, optimizer, callback=callback, maxiters=max_iterations,
                abstol=settings.inversion_settings.inversion_abs_tols[iOptimizer], reltol=settings.inversion_settings.inversion_rel_tols[iOptimizer])
        elseif optimizer_choice == "LBFGS"   #LBFGS does not support abstol
            sol = solve(optprob, optimizer, callback=callback, maxiters=max_iterations, reltol=settings.inversion_settings.inversion_rel_tols[iOptimizer])
        else
            error("Not implemented yet")
        end

        #check whether the optimizer converged. If converged, break the loop (no need to continue)
        if sol.retcode == SciMLBase.ReturnCode.Success
            println("   Optimizer $(optimizer_choice) converged. No need to continue.")
            break
        end

        #Manually do the training loop: Have more control on the optimizer
        #for epoch in 1:settings.UDE_settings.UDE_max_iterations
        # Compute gradients and update parameters
        #    grads = Zygote.pullback(θ -> compute_loss_UDE(ode_f, Q0, tspan, θ, settings, my_mesh_2D, swe_2D_constants, observed_data, data_type), θ)
        #    θ = optimizer(θ, grads)
        #end

        #update the initial condition for the next optimizer
        optprob = OptimizationProblem(optf, copy(sol.u))  # No parameters needed
        GC.gc() #clear the memory

    end

    return sol, ITER, LOSS, PRED, PARS
end

# Bound loss function compatible with both Zygote and ForwardDiff
function compute_bound_loss(params::Vector{T}, lower_bound::Real, upper_bound::Real) where {T}
    loss = zero(T)
    lb = convert(float(T), lower_bound)  # Convert bounds to parameter type
    ub = convert(float(T), upper_bound)
    for p in params
        loss += sum(abs2, max(zero(T), lb - p)) +  # Lower bound violation
                sum(abs2, max(zero(T), p - ub))    # Upper bound violation
    end
    return loss
end

#loss due to slope (gradient) regularization
function calc_slope_loss(params::Vector{T}, settings, my_mesh_2D) where {T}

    #make sure params is a field (defined at cell centers)
    if length(params) != my_mesh_2D.numOfCells
        throw(ArgumentError("params must be a field (defined at cell centers)"))
    end

    # Compute gradient at cell centers 
    S0 = compute_scalar_gradients(my_mesh_2D, params)

    #compute the magnitude of S0
    S0_mag = sqrt.(S0[:, 1] .^ 2 + S0[:, 2] .^ 2 .+ eps(T))

    #compute the loss due to slope regularization
    loss = compute_bound_loss(S0_mag, 0.0, settings.inversion_settings.parameter_slope_upper_bound)

    return loss
end