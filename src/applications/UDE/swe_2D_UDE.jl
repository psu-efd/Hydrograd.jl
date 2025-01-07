
#perform UDE for 2D SWE (for model calibration/parameter estimation/hidden physics discovery/etc.)
using Dates

#using Debugger, JuliaInterpreter

#main function for the UDE
function swe_2D_UDE(ode_f, Q0, params_vector, swe_extra_params)

    #unpack the extra parameters
    settings = swe_extra_params.settings
    my_mesh_2D = swe_extra_params.my_mesh_2D
    swe_2D_constants = swe_extra_params.swe_2D_constants
    case_path = swe_extra_params.case_path

    nodeCoordinates = swe_extra_params.nodeCoordinates

    if settings.bVerbose
        println("UDE choice: ", settings.UDE_settings.UDE_choice)
    end

    #open the forward simulation result (as the ground truth)
    sol_truth = JSON3.read(open(joinpath(case_path, settings.UDE_settings.UDE_truth_file_name)), Dict{String,Vector{Float64}})

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
    ManningN_cell_truth = vec(sol_truth["ManningN_cell_truth"])
    inlet_discharges_truth = vec(sol_truth["inlet_discharges_truth"])

    println("   Loading calibration data ...")

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
        end
    end

    #p_init is the initial weights for the UDE's NN model
    p_init = ComponentVector{Float64}(swe_extra_params.ude_model_params)

    #@show typeof(p_init)
    #@show p_init

    #perform the UDE
    sol, ITER, LOSS, PRED, PARS = nothing, nothing, nothing, nothing, nothing
    if settings.UDE_settings.UDE_mode == "training"
        println("   Performing UDE training ...\n")
        sol, ITER, LOSS, PRED, PARS = UDE_training(ode_f, Q0, swe_2D_constants.tspan, p_init, settings, my_mesh_2D, swe_2D_constants, observed_data, case_path)

        #@show typeof(sol.u)
        #@show sol.u

        #save the UDE results
        jldsave(joinpath(case_path, settings.UDE_settings.UDE_training_save_file_name); ITER, LOSS, PRED, PARS,
            ude_model_params=sol.u, ude_model_state=swe_extra_params.ude_model_state)

        #process the UDE results
        println("   Post-processing UDE results ...")

        #process UDE results
        #so.u is the trained NN weights (ude_model_params)
        Hydrograd.postprocess_UDE_training_results_swe_2D(swe_extra_params, zb_cell_truth, h_truth, u_truth, v_truth, WSE_truth)
    elseif settings.UDE_settings.UDE_mode == "inference"
        println("   Performing UDE inference ...\n")
        sol = UDE_inference(ode_f, Q0, swe_2D_constants.tspan, p_init, settings, my_mesh_2D, swe_2D_constants, observed_data, case_path)

        #save the UDE results
        jldsave(joinpath(case_path, settings.UDE_settings.UDE_inference_save_file_name); sol)

        #process the UDE results
        println("   Post-processing UDE results ...")

        #process UDE results
        Hydrograd.postprocess_UDE_inference_results_swe_2D(settings, my_mesh_2D, nodeCoordinates, zb_cell_truth, h_truth, u_truth, v_truth, WSE_truth, case_path)
    else
        error("Invalid UDE mode: $(settings.UDE_settings.UDE_mode). Supported options: training, inference.")
    end

end

#train the UDE
function UDE_training(ode_f, Q0, tspan, p_init, settings, my_mesh_2D, swe_2D_constants, observed_data, case_path)
    #start the timer for the UDE
    UDE_start_time = now()  # Current date and time

    # Loss function for optimization: θ is the trainable parameters (NN weights), p is not used
    function opt_loss(θ, p)  # Add p argument even if unused

        data_type = eltype(θ)

        Zygote.ignore() do
            #println("data_type = ", data_type)
            #println("θ = ", θ)
            #println("p_init = ", p_init)
        end

        loss_total, loss_pred_WSE, loss_pred_uv, pred = compute_loss_UDE(ode_f, Q0, tspan, θ, settings,
            my_mesh_2D, swe_2D_constants, observed_data, data_type)

        # Call callback with all values (but outside gradient calculation)
        Zygote.ignore() do
            if settings.UDE_settings.UDE_optimizer == "LBFGS"
                callback_within_opt_loss(θ, loss_total, loss_pred_WSE, loss_pred_uv, pred)
            end
        end

        #return loss_total, loss_pred_WSE, loss_pred_uv, pred
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
    if settings.UDE_settings.UDE_sensealg == "AutoZygote()"
        adtype = Optimization.AutoZygote()
    elseif settings.UDE_settings.UDE_sensealg == "AutoReverseDiff()"
        adtype = Optimization.AutoReverseDiff(compile=false)
    elseif settings.UDE_settings.UDE_sensealg == "AutoForwardDiff()"
        adtype = Optimization.AutoForwardDiff()
    elseif settings.UDE_settings.UDE_sensealg == "AutoReverseDiff()"
        adtype = Optimization.AutoReverseDiff(compile=false)
    else
        println("       settings.UDE_settings.UDE_sensealg = ", settings.UDE_settings.UDE_sensealg)
        throw(ArgumentError("Invalid sensealg choice. Supported sensealg: AutoZygote(), AutoReverseDiff(), AutoForwardDiff(), AutoFiniteDiff(), AutoModelingToolkit(), AutoEnzyme(). No UDE is performed."))
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

    #create the optimizer
    if settings.UDE_settings.UDE_optimizer == "Adam"
        optimizer = Adam(settings.UDE_settings.UDE_learning_rate)
    elseif settings.UDE_settings.UDE_optimizer == "LBFGS"
        optimizer = Optimization.LBFGS()
    else
        throw(ArgumentError("Invalid optimizer choice. Supported optimizers: Adam, LBFGS. No UDE is performed."))
    end

    #define the accumulators for the UDE results
    ITER = []        #iteration number accumulator
    LOSS = []        # Loss accumulator
    PRED = []        # prediction accumulator
    PARS = []        # parameters accumulator

    # Define the callback 
    callback_within_opt_loss = function (θ, loss_total, loss_pred_WSE, loss_pred_uv, prediction) #callback function to observe training
        UDE_current_time = now()  # Current date and time
        UDE_elapsed_time = UDE_current_time - UDE_start_time
        UDE_elapsed_seconds = Millisecond(UDE_elapsed_time).value / 1000
        UDE_start_time = UDE_current_time

       iter_number = length(ITER) + 1

       append!(ITER, iter_number)
       append!(LOSS, [[loss_total, loss_pred_WSE, loss_pred_uv]])
       append!(PRED, [prediction[:, 1, end]])
       append!(PARS, [θ])
       
       println("iter_number = ", iter_number, ", loss_total = ", loss_total, ", UDE_elapsed_seconds = ", UDE_elapsed_seconds)

       return false
   end

    callback = function (optimizer_state, loss_total) #callback function to observe training
         UDE_current_time = now()  # Current date and time
         UDE_elapsed_time = UDE_current_time - UDE_start_time
         UDE_elapsed_seconds = Millisecond(UDE_elapsed_time).value / 1000
         UDE_start_time = UDE_current_time

        θ = optimizer_state.u

        iter_number = length(ITER) + 1

        append!(ITER, iter_number)
        append!(LOSS, [[loss_total]])

        append!(PARS, [θ])
        
        println("iter_number = ", iter_number, ", loss_total = ", loss_total, ", UDE_elapsed_seconds = ", UDE_elapsed_seconds)

        return false
    end
    # The first argument is the state of the optimizer, which is an OptimizationState object
    # callback = function (optimizer_state, loss_total, loss_pred_WSE, loss_pred_uv, prediction) #callback function to observe training

    #     # Extract the parameters from the state
    #     θ = optimizer_state.u

    #     #@show typeof(θ)
    #     #@show θ

    #     # Print gradient norms per layer
    #     # gs = first(Zygote.gradient(θ) do p
    #     #     loss_total, _, _, _ = opt_loss(p, nothing)
    #     #     return loss_total
    #     # end)

    #     # println("Structure of gs: ", keys(gs))
    #     # for (k, v) in pairs(gs)
    #     #     println("Key: ", k)
    #     #     println("Value type: ", typeof(v))
    #     #     if v isa NamedTuple
    #     #         println("  Subfields: ", keys(v))
    #     #     end
    #     # end
        
    #     # # Analyze gradients by layer
    #     # println("\nGradient analysis by layer:")
        
    #     # # Layer 1 gradients
    #     # layer1_weights = gs.layer_1.weight
    #     # layer1_bias = gs.layer_1.bias
    #     # println("Layer 1:")
    #     # println("  weights gradient = ", layer1_weights)
    #     # println("  bias gradient = ", layer1_bias)
    #     # println("  Weights gradient norm: ", sqrt(sum(abs2, layer1_weights)))
    #     # println("  Bias gradient norm: ", sqrt(sum(abs2, layer1_bias)))
        
    #     # # Layer 2 gradients
    #     # layer3_weights = gs.layer_3.weight
    #     # layer3_bias = gs.layer_3.bias
    #     # println("Layer 2:")
    #     # println("  weights gradient = ", layer3_weights)
    #     # println("  bias gradient = ", layer3_bias)
    #     # println("  Weights gradient norm: ", sqrt(sum(abs2, layer3_weights)))
    #     # println("  Bias gradient norm: ", sqrt(sum(abs2, layer3_bias)))
        
    #     # # Layer 3 gradients
    #     # layer5_weights = gs.layer_5.weight
    #     # layer5_bias = gs.layer_5.bias
    #     # println("Layer 3:")
    #     # println("  weights gradient = ", layer5_weights)
    #     # println("  bias gradient = ", layer5_bias)
    #     # println("  Weights gradient norm: ", sqrt(sum(abs2, layer5_weights)))
    #     # println("  Bias gradient norm: ", sqrt(sum(abs2, layer5_bias)))

    #     UDE_current_time = now()  # Current date and time
    #     UDE_elapsed_time = UDE_current_time - UDE_start_time
    #     UDE_elapsed_seconds = Millisecond(UDE_elapsed_time).value / 1000

    #     #reset the UDE start time
    #     UDE_start_time = UDE_current_time

    #     iter_number = size(LOSS)[1]  #get the UDE iteration number (=length of LOSS array)

    #     Zygote.ignore() do
    #         println("      iter, loss_total, loss_pred_WSE, loss_pred_uv = ", iter_number, ", ",
    #             loss_total, ", ", loss_pred_WSE, ", ", loss_pred_uv)
    #         println("      UDE iteration elapsed seconds = ", UDE_elapsed_seconds)
    #     end

    #     #save the UDE results
    #     if iter_number % settings.UDE_settings.UDE_training_save_frequency == 0
    #         append!(ITER, iter_number)

    #         append!(PRED, [prediction[:, 1, end]])

    #         append!(LOSS, [[loss_total, loss_pred_WSE, loss_pred_uv]])

    #         append!(PARS, [θ])
    #     end

    #     #checkpoint the UDE results (in case the UDE training is interrupted)
    #     if settings.UDE_settings.UDE_training_save_checkpoint &&
    #        iter_number % settings.UDE_settings.UDE_training_checkpoint_frequency == 0 &&
    #        iter_number > 0
    #         jldsave(joinpath(case_path, "checkpoint_UDE_iter_$(iter_number).jld2"); ITER, LOSS, PRED, PARS)
    #     end

    #     #if l > 1e-9
    #     #    false
    #     #else 
    #     #    true   #force the optimizer to stop 
    #     #end

    #     false
    # end


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

    #Solve in one call
    if settings.UDE_settings.UDE_optimizer == "Adam"
        sol = solve(optprob, optimizer, callback=callback, maxiters=settings.UDE_settings.UDE_max_iterations, abstol=settings.UDE_settings.UDE_abs_tol, reltol=settings.UDE_settings.UDE_rel_tol)
    elseif settings.UDE_settings.UDE_optimizer == "LBFGS"
        sol = solve(optprob, optimizer, maxiters=settings.UDE_settings.UDE_max_iterations, reltol=settings.UDE_settings.UDE_rel_tol)
    else
        error("UDE optimizer not implemented")
    end

    #Manually do the training loop: Have more control on the optimizer
    #for epoch in 1:settings.UDE_settings.UDE_max_iterations
        # Compute gradients and update parameters
    #    grads = Zygote.pullback(θ -> compute_loss_UDE(ode_f, Q0, tspan, θ, settings, my_mesh_2D, swe_2D_constants, observed_data, data_type), θ)
    #    θ = optimizer(θ, grads)
    #end
    

    return sol, ITER, LOSS, PRED, PARS
end


# Define the loss function
function compute_loss_UDE(ode_f, Q0, tspan, p, settings, my_mesh_2D, swe_2D_constants, observed_data, data_type)
    # Solve the ODE (forward pass)
    pred = UDE_forward_simulation(ode_f, Q0, tspan, p, settings, swe_2D_constants)

    #compute the loss
    # Ensure type stability in loss computation
    loss_total = zero(data_type)
    loss_pred_WSE = zero(data_type)
    loss_pred_uv = zero(data_type)

    if pred.retcode == SciMLBase.ReturnCode.Success
        WSE_truth = Vector{Float64}(observed_data["WSE_truth"])
        #h_truth = Vector{Float64}(observed_data["h_truth"])
        u_truth = Vector{Float64}(observed_data["u_truth"])
        v_truth = Vector{Float64}(observed_data["v_truth"])

        zb_cells_temp = observed_data["zb_cell_truth"]

        l = pred[:, 1, end] .+ zb_cells_temp .- WSE_truth  #loss = free surface elevation mismatch

        #loss for free surface elevation mismatch
        if settings.UDE_settings.UDE_bWSE_loss
            loss_pred_WSE = sum(abs2, l)
        end

        #loss for velocity mismatch
        # Add small epsilon to prevent division by zero
        ϵ = sqrt(eps(data_type))

        if settings.UDE_settings.UDE_b_uv_loss      #if also include u in the loss 
            l_u = pred[:, 2, end] ./ (pred[:, 1, end] .+ ϵ) .- u_truth
            l_v = pred[:, 3, end] ./ (pred[:, 1, end] .+ ϵ) .- v_truth

            loss_pred_uv = sum(abs2, l_u) + sum(abs2, l_v)
        end

        #combined loss due to free surface elevation mismatch and velocity mismatch
        loss_total = loss_pred_WSE + loss_pred_uv
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

    return loss_total, loss_pred_WSE, loss_pred_uv, pred
end



#inference from the trained UDE
function UDE_inference(ode_f, Q0, tspan, p_init, settings, my_mesh_2D, swe_2D_constants, observed_data, case_path)

end

#perform forward simulation
function UDE_forward_simulation(ode_f, Q0, tspan, p, settings, swe_2D_constants)

    # Create ODEProblem: p should be the weights of the NN model         
    prob = ODEProblem(ode_f, Q0, tspan, p)

    # Create ODE solver
    if settings.UDE_settings.UDE_ode_solver == "Tsit5()"
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
    if settings.UDE_settings.UDE_ode_solver_sensealg == "ZygoteVJP()"
        ode_solver_sensealg = ZygoteVJP()
    elseif settings.UDE_settings.UDE_ode_solver_sensealg == "InterpolatingAdjoint(autojacvec=ZygoteVJP())"
        ode_solver_sensealg = InterpolatingAdjoint(autojacvec=ZygoteVJP())
    elseif settings.UDE_settings.UDE_ode_solver_sensealg == "ForwardDiffSensitivity()"
        ode_solver_sensealg = ForwardDiffSensitivity()
    elseif settings.UDE_settings.UDE_ode_solver_sensealg == "BacksolveAdjoint(autojacvec=ZygoteVJP())"
        ode_solver_sensealg = BacksolveAdjoint(autojacvec=ZygoteVJP())
    else
        error("Not implemented yet")
    end

    #define the time for saving the results for the ODE solver
    dt_save = (swe_2D_constants.tspan[2] - swe_2D_constants.tspan[1]) / settings.UDE_settings.UDE_ode_solver_nSave
    t_save = swe_2D_constants.tspan[1]:dt_save:swe_2D_constants.tspan[2]

    # Solve the ODE
    pred = solve(prob, ode_solver, adaptive=settings.UDE_settings.UDE_ode_solver_adaptive, dt=swe_2D_constants.dt,
        saveat=t_save, sensealg=ode_solver_sensealg)

    Zygote.ignore() do
        if settings.bVerbose
            #@show pred.retcode
            #@show typeof(pred)
            #@show size(pred)
            #@show pred
        end
    end

    return pred
end
