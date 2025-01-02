
#perform inversion for 2D SWE

#using Debugger, JuliaInterpreter

# Define the loss function
function compute_loss(ode_f, Q0, tspan, p, settings, my_mesh_2D, swe_2D_constants, observed_data, active_param_name, data_type)
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
    elseif settings.inversion_settings.ode_solver_sensealg == "ForwardDiffSensitivity()"
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
        loss_pred_WSE = sum(abs2, l)

        #loss for velocity mismatch
        # Add small epsilon to prevent division by zero
        ϵ = sqrt(eps(data_type))

        if settings.inversion_settings.bInversion_u_loss      #if also include u in the loss 
            l_u = pred[:, 2, end] ./ (pred[:, 1, end] .+ ϵ) .- u_truth
            l_v = pred[:, 3, end] ./ (pred[:, 1, end] .+ ϵ) .- v_truth

            loss_pred_uv = sum(abs2, l_u) + sum(abs2, l_v)
        end

        #combined loss due to free surface elevation mismatch and velocity mismatch
        loss_pred = loss_pred_WSE + loss_pred_uv

        #loss for bed slope regularization
        if settings.inversion_settings.bInversion_slope_loss    #if bed slope is included in the loss 
            loss_slope = calc_slope_loss(zb_cells_temp, my_mesh_2D)
        end

        #combined loss due to free surface elevation mismatch, velocity mismatch, and bed slope regularization
        loss_total = loss_pred + loss_slope
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

    return loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope, pred
end

function optimize_parameters(ode_f, Q0, tspan, p_init, settings, my_mesh_2D, swe_2D_constants, observed_data, active_param_name)
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

        loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope, pred = compute_loss(ode_f, Q0, tspan, θ, settings, 
                                                                                            my_mesh_2D, swe_2D_constants, observed_data, active_param_name, data_type)

        # Call callback with all values (but outside gradient calculation)
        Zygote.ignore() do
            #callback(θ, loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope, pred)
        end

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

    #create the optimizer
    if settings.inversion_settings.optimizer == "Adam"
        optimizer = Adam(settings.inversion_settings.learning_rate)
    elseif settings.inversion_settings.optimizer == "LBFGS"
        optimizer = LBFGS()
    else
        throw(ArgumentError("Invalid optimizer choice. Supported optimizers: Adam, LBFGS. No inversion is performed."))
    end

     # Define the callback to handle both vector and OptimizationState inputs
     callback = function (state_or_θ, loss)
        θ = state_or_θ isa Optimization.OptimizationState ? state_or_θ.u : state_or_θ
        println("Loss: ", loss)
        println("Parameters: ", θ)
        return false  # continue optimization
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

    sol = solve(optprob, optimizer, callback=callback, maxiters=settings.inversion_settings.max_iterations)
    #sol = solve(optprob, Adam(settings.inversion_settings.learning_rate), maxiters=settings.inversion_settings.max_iterations)
    #@enter sol = @interpret solve(optprob, Adam(settings.inversion_settings.learning_rate), callback=callback, maxiters=settings.inversion_settings.max_iterations)

    return sol
end

#main function for the inversion
function swe_2D_inversion(ode_f, Q0, params_vector, swe_extra_params, case_path)

    #unpack the extra parameters
    active_param_name = swe_extra_params.active_param_name
    settings = swe_extra_params.settings
    my_mesh_2D = swe_extra_params.my_mesh_2D
    swe_2D_constants = swe_extra_params.swe_2D_constants

    #if settings.bVerbose
        println("       swe_2D_inversion")
        println("       inversion parameter name = ", active_param_name)
    #end

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
    ManningN_cell_truth = vec(sol_truth["ManningN_cell_truth"])
    inlet_discharges_truth = vec(sol_truth["inlet_discharges_truth"])

    if settings.bVerbose
        #@show WSE_truth
        #@show h_truth
        #@show u_truth
        #@show v_truth
        #@show zb_cell_truth
        #@show ManningN_cell_truth
        #@show inlet_discharges_truth
    end

    println("       Loading inversion data ...")

    #combine the truth data into a dictionary
    observed_data = Dict{String,Vector{Float64}}("WSE_truth" => WSE_truth, "h_truth" => h_truth, "u_truth" => u_truth, "v_truth" => v_truth,
        "zb_cell_truth" => zb_cell_truth, "ManningN_cell_truth" => ManningN_cell_truth, "inlet_discharges_truth" => inlet_discharges_truth)

    #debug AD correctness
    #debug start
    #debug_AD(ode_f, Q0, swe_2D_constants, params_vector, swe_extra_params)
    #return
    #debug end


    #define the accumulators for the inversion results
    LOSS::Vector{Float64} = []                              # Loss accumulator
    #PRED::Vector{Any} = []                              # prediction accumulator
    PARS::Vector{Vector{Float64}} = []                              # parameters accumulator

    
   
    # callback = let
    #     PARS = PARS
    #     LOSS = LOSS
    #     (θ::Vector{Float64}, loss_total::Float64) -> begin
    #         iter = length(LOSS)
    #         Zygote.ignore() do
    #             println("      iter, loss_total = ", iter, ", ", loss_total)
    #         end
    #         push!(LOSS, loss_total)
    #         push!(PARS, copy(θ))
    #         false
    #     end
    # end

    # callback = function (θ, loss_total) #callback function to observe training
    #         iter = size(LOSS)[1]  #get the inversion iteration number (=length of LOSS array)

    #         Zygote.ignore() do
    #             println("      iter, loss_total = ", iter, ", ", loss_total)
    #         end

    #         #append!(PRED, [pred[:, 1, end]])
    #         append!(LOSS, [[loss_total]])

    #         if !isa(θ, Vector{Float64})  #NLopt returns an optimization object, not an arrary
    #             #println("theta.u = ", θ.u)
    #             append!(PARS, [copy(θ.u)])
    #         else
    #             append!(PARS, θ)
    #         end

    #         #if l > 1e-9
    #         #    false
    #         #else 
    #         #    true   #force the optimizer to stop 
    #         #end

    #         false
    #     end

    # callback = function (θ, loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope, pred) #callback function to observe training
    #     iter = size(LOSS)[1]  #get the inversion iteration number (=length of LOSS array)

    #     Zygote.ignore() do
    #         println("      iter, loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope = ", iter, ", ",
    #             loss_total, ", ", loss_pred, ", ", loss_pred_WSE, ", ", loss_pred_uv, ", ", loss_slope)
    #     end

    #     append!(PRED, [pred[:, 1, end]])
    #     append!(LOSS, [[loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_slope, pred]])

    #     if !isa(θ, Vector{Float64})  #NLopt returns an optimization object, not an arrary
    #         #println("theta.u = ", θ.u)
    #         append!(PARS, [copy(θ.u)])
    #     else
    #         append!(PARS, θ)
    #     end

    #     #if l > 1e-9
    #     #    false
    #     #else 
    #     #    true   #force the optimizer to stop 
    #     #end

    #     false
    # end

    Zygote.ignore() do
        @show typeof(Q0)
        @show typeof(swe_2D_constants.tspan)
        @show typeof(observed_data)
        @show typeof(params_vector)

        @show size(params_vector)
        @show params_vector
    end

    #perform the inversion
    println("       Performing inversion ...")
    sol = optimize_parameters(ode_f, Q0, swe_2D_constants.tspan, params_vector, settings, my_mesh_2D, swe_2D_constants, observed_data, active_param_name)

    #save the inversion results
    #jldsave(joinpath(case_path, settings.inversion_settings.save_file_name); LOSS, PRED, PARS)

    #process the inversion results
    #println("   Post-processing inversion results ...")

    #process inversion results
    #Hydrograd.postprocess_inversion_results_swe_2D(settings, my_mesh_2D, nodeCoordinates, zb_cell_truth, h_truth, u_truth, v_truth, WSE_truth, case_path)

end
