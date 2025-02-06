
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

    #Note: ManningN_cells is updated in the semi_discretize_swe_2D function. Thus, ManningN_cells in swe2d_extra_params cannot be used if the UDE is for n(h) (because it can not be updated in the imutable struct swe2d_extra_params).
    ManningN_cells = swe_extra_params.ManningN_cells

    wstill = swe_extra_params.wstill
    hstill = swe_extra_params.hstill

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
    ManningN_cells_truth = vec(sol_truth["ManningN_cells_truth"])
    friction_x_truth = vec(sol_truth["friction_x_truth"])
    friction_y_truth = vec(sol_truth["friction_y_truth"])
    inlet_discharges_truth = vec(sol_truth["inlet_discharges_truth"])

    println("   Loading calibration data ...")

    #combine the truth data into a dictionary
    observed_data = Dict{String,Vector{Float64}}("WSE_truth" => WSE_truth, "h_truth" => h_truth, "u_truth" => u_truth, "v_truth" => v_truth,
        "zb_cell_truth" => zb_cell_truth, "ManningN_cells_truth" => ManningN_cells_truth,
        "friction_x_truth" => friction_x_truth, "friction_y_truth" => friction_y_truth,
        "inlet_discharges_truth" => inlet_discharges_truth)

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
    sol, ITER, LOSS, PARS = nothing, nothing, nothing, nothing
    if settings.UDE_settings.UDE_mode == "training"
        println("   Performing UDE training: ", settings.UDE_settings.UDE_choice, " ...\n")
        sol, ITER, LOSS, PARS = UDE_training(ode_f, Q0, swe_2D_constants.tspan, p_init, settings, my_mesh_2D, swe_2D_constants, observed_data, wstill, hstill, case_path)

        #@show typeof(sol.u)
        #@show sol.u

        #save the UDE results
        jldsave(joinpath(case_path, settings.UDE_settings.UDE_training_save_file_name); ITER, LOSS, PARS,
            ude_model_params=sol.u, ude_model_state=swe_extra_params.ude_model_state)

        #compute all the losses and ODE solutions for each inversion iteration
        println("   Computing all the losses and ODE solutions for each UDE iteration ...")
        compute_all_losses_UDE(ode_f, Q0, settings, swe_extra_params, nodeCoordinates, my_mesh_2D, swe_2D_constants, observed_data, wstill, hstill, ManningN_cells, case_path)        

        #process the UDE results
        println("   Post-processing UDE results ...")

        #process UDE results
        #so.u is the trained NN weights (ude_model_params)
        #Hydrograd.postprocess_UDE_training_results_swe_2D(swe_extra_params, zb_cell_truth, h_truth, u_truth, v_truth, WSE_truth,
        #    ManningN_cells_truth, friction_x_truth, friction_y_truth)

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
function UDE_training(ode_f, Q0, tspan, p_init, settings, my_mesh_2D, swe_2D_constants, observed_data, wstill, hstill, case_path)
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

        loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, pred = compute_loss_UDE(ode_f, Q0, tspan, θ, settings,
            my_mesh_2D, swe_2D_constants, observed_data, wstill, hstill, data_type)

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

    #define the accumulators for the UDE results
    ITER = []        #iteration number accumulator
    LOSS = []        # Loss accumulator
    PARS = []        # parameters accumulator

    # Define the callback 
    callback = function (optimizer_state, loss_total) #callback function to observe training
        UDE_current_time = now()  # Current date and time
        UDE_elapsed_time = UDE_current_time - UDE_start_time
        UDE_elapsed_seconds = Millisecond(UDE_elapsed_time).value / 1000

        #reset the UDE start time
        UDE_start_time = UDE_current_time

        iter_number = length(ITER) + 1

        θ = optimizer_state.u

        #append the loss and parameters to the accumulators
        append!(ITER, iter_number)
        append!(LOSS, [loss_total])
        append!(PARS, [θ])

        Zygote.ignore() do
            println("       iter_number = ", iter_number, ", loss_total = ", loss_total, ", UDE_elapsed_seconds = ", UDE_elapsed_seconds)

            #if not the first iteration, compute and print the max and min of the parameter change from the previous iteration
            if iter_number > 1
                param_change = θ - PARS[end-1]
                println("           max(param_change) = ", maximum(param_change), ", min(param_change) = ", minimum(param_change))
            end

            #save the UDE results (both JSON and JLD2; a little bit redundant but just for the sake of safety)
            if iter_number % settings.UDE_settings.UDE_training_save_frequency == 0

                data_type = eltype(θ)
                
                call_back_save_dict = Dict("iter_number" => iter_number, "loss_total" => loss_total, "theta" => θ)
                open(joinpath(case_path, "UDE_callback_save_iter_$(iter_number).json"), "w") do io
                    JSON3.pretty(io, call_back_save_dict)
                end

                jldsave(joinpath(case_path, "UDE_callback_save_iter_$(iter_number).jld2"); iter_number=iter_number, loss_total=loss_total, theta=θ)
            end

            #checkpoint the UDE results (in case the UDE training is interrupted)
            if settings.UDE_settings.UDE_training_save_checkpoint &&
               iter_number % settings.UDE_settings.UDE_training_checkpoint_frequency == 0 &&
               iter_number > 0
                jldsave(joinpath(case_path, "checkpoint_UDE_iter_$(iter_number).jld2"); ITER, LOSS, PARS)
            end

        end

        return false
    end

    sol = nothing

    #loop over the optimizers
    for (iOptimizer, optimizer_choice) in enumerate(settings.UDE_settings.UDE_optimizers)
        println("   iOptimizer: ", iOptimizer, " out of ", length(settings.UDE_settings.UDE_optimizers))
        println("       optimizer_choice: ", optimizer_choice)
        println("       iterations: ", settings.UDE_settings.UDE_max_iterations[iOptimizer])
        println("       learning rate: ", settings.UDE_settings.UDE_learning_rates[iOptimizer])
        println("       abs_tol: ", settings.UDE_settings.UDE_abs_tols[iOptimizer])
        println("       rel_tol: ", settings.UDE_settings.UDE_rel_tols[iOptimizer])

        max_iterations = settings.UDE_settings.UDE_max_iterations[iOptimizer]

        if max_iterations <= 0
            println("       max_iterations <= 0. No UDE is performed with this optimizer.")
            continue
        end

        #create the optimizer
        if optimizer_choice == "Adam"
            optimizer = OptimizationOptimisers.Adam(settings.UDE_settings.UDE_learning_rates[iOptimizer])
        elseif optimizer_choice == "LBFGS"
            optimizer = OptimizationOptimJL.LBFGS(linesearch = LineSearches.BackTracking())
        else
            throw(ArgumentError("Invalid optimizer choice. Supported optimizers: Adam, LBFGS. No UDE is performed."))
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

        #Solve in one call
        if optimizer_choice == "Adam"
            sol = solve(optprob, optimizer, callback=callback, maxiters=settings.UDE_settings.UDE_max_iterations[iOptimizer], 
                         abstol=settings.UDE_settings.UDE_abs_tols[iOptimizer], reltol=settings.UDE_settings.UDE_rel_tols[iOptimizer])
        elseif optimizer_choice == "LBFGS"  #LBFGS does not support abstol
            sol = solve(optprob, optimizer, callback=callback, maxiters=settings.UDE_settings.UDE_max_iterations[iOptimizer], reltol=settings.UDE_settings.UDE_rel_tols[iOptimizer])
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


    return sol, ITER, LOSS, PARS
end

# compute all the losses for each UDE training iteration
function compute_all_losses_UDE(ode_f::ODEFunction, Q0::AbstractVector{T1}, settings::ControlSettings, swe_extra_params::SWE2D_Extra_Parameters, nodeCoordinates::Matrix{T2}, my_mesh_2D::mesh_2D, swe_2D_constants::swe_2D_consts, observed_data::Dict{String,Vector{T3}}, 
    wstill::AbstractVector{T4}, hstill::AbstractVector{T5}, ManningN_cells::AbstractVector{T6}, case_path::String) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real,T6<:Real}

    #Note: ManningN_cells needs to be updated from h before use. 

    data_type = eltype(Q0)
    
    #find the list of files of "UDE_callback_save_iter_<iter_number>.jld2"
    iter_files = filter(file -> occursin("UDE_callback_save_iter_", file) && endswith(file, ".jld2"), readdir(case_path))

    #println("iter_files = ", iter_files)

    #sort the files by the iteration number
    iter_files = sort(iter_files, by = file -> begin
        # Extract the number between "iter_" and ".json"
        m = match(r"UDE_callback_save_iter_(\d+)\.jld2$", file)
        if isnothing(m)
            error("Invalid filename format: $file")
        end
        parse(Int, m.captures[1])
    end)

    #print the total number of iterations
    nIter = length(iter_files)
    println("       Total number of UDE iterations: ", nIter)

    #get the truth data
    WSE_truth = observed_data["WSE_truth"]
    h_truth = observed_data["h_truth"]
    u_truth = observed_data["u_truth"]
    v_truth = observed_data["v_truth"]
    zb_cell_truth = observed_data["zb_cell_truth"]
    ManningN_cell_truth = observed_data["ManningN_cells_truth"]
    
    #initialize the losses
    loss_totals = []
    loss_preds = []
    loss_pred_WSEs = []
    loss_pred_uvs = []
    inverted_WSEs = []
    inverted_hs = []
    inverted_us = []
    inverted_vs = []
    
    #load the files
    for file in iter_files
        println("       Loading UDE results from file: ", file)
        data = load(joinpath(case_path, file))     

        iter_number = parse(Int, match(r"UDE_callback_save_iter_(\d+)\.jld2$", file).captures[1])
        iter_number_from_file = data["iter_number"]
        @assert iter_number == iter_number_from_file "The iteration number in the file name and the file content do not match: iter_number = $iter_number, iter_number_from_file = $iter_number_from_file"

        θ_from_file = data["theta"]
        total_loss_from_file = convert(Float64, data["loss_total"])

        #@show typeof(θ_from_file)

        loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, ode_pred = compute_loss_UDE(ode_f, Q0, swe_2D_constants.tspan, θ_from_file, settings, 
            my_mesh_2D, swe_2D_constants, observed_data, wstill, hstill, data_type)

        #assert the total loss is the same as the one in the file
        @assert abs(loss_total - total_loss_from_file) < 1e-3 "The total loss in the file and the computed loss do not match: loss_total = $loss_total, total_loss_from_file = $total_loss_from_file"

        #append the losses to the accumulators
        push!(loss_totals, loss_total)
        push!(loss_preds, loss_pred)
        push!(loss_pred_WSEs, loss_pred_WSE)
        push!(loss_pred_uvs, loss_pred_uv)

        zb_i = zb_cell_truth

        #@show size(ode_pred.u[end])
        #@show size(zb_i)

        #append the ODE solution to the accumulator
        push!(inverted_WSEs, ode_pred.u[end][1:my_mesh_2D.numOfCells] .+ hstill .+ zb_i)
        push!(inverted_hs, ode_pred.u[end][1:my_mesh_2D.numOfCells] .+ hstill)
        push!(inverted_us, ode_pred.u[end][my_mesh_2D.numOfCells+1:2*my_mesh_2D.numOfCells] ./ (ode_pred.u[end][1:my_mesh_2D.numOfCells] .+ hstill)) 
        push!(inverted_vs, ode_pred.u[end][2*my_mesh_2D.numOfCells+1:3*my_mesh_2D.numOfCells] ./ (ode_pred.u[end][1:my_mesh_2D.numOfCells] .+ hstill))

        #save the inverted results to vtk
        bSave_vtk = true
        if bSave_vtk
            println("       Saving UDE iteration results to vtk file for iteration: ", iter_number)

            field_name = "iter_number"
            field_type = "integer"
            field_value = iter_number

            xi_array = ode_pred.u[end][1:my_mesh_2D.numOfCells]
            h_array = xi_array .+ hstill
            q_x_array = ode_pred.u[end][my_mesh_2D.numOfCells+1:2*my_mesh_2D.numOfCells]
            q_y_array = ode_pred.u[end][2*my_mesh_2D.numOfCells+1:3*my_mesh_2D.numOfCells]

            u_temp = q_x_array ./ (h_array .+ eps(eltype(h_array)))
            v_temp = q_y_array ./ (h_array .+ eps(eltype(h_array)))
            U_vector = hcat(u_temp, v_temp)
            Umag = sqrt.(u_temp.^2 .+ v_temp.^2)
                            
            vector_data = [U_vector] 
            vector_names = ["U"]

            WSE_array = h_array + zb_i

            #update ManningN_cells from h and the NN if needed
            if settings.UDE_settings.UDE_choice == "ManningN_h" || settings.UDE_settings.UDE_choice == "ManningN_h_Umag_ks"
                ManningN_cells = update_ManningN_UDE(settings.UDE_settings.UDE_choice, 
                           h_array, Umag, swe_extra_params.ks_cells, 
                           swe_extra_params.ude_model, θ_from_file, swe_extra_params.ude_model_state, 
                           Float64.(settings.UDE_settings.UDE_NN_config["h_bounds"]), 
                           Float64.(settings.UDE_settings.UDE_NN_config["Umag_bounds"]), 
                           Float64.(settings.UDE_settings.UDE_NN_config["ks_bounds"]), 
                           my_mesh_2D.numOfCells)                        
            end

            ks_cells = swe_extra_params.ks_cells
            h_ks_cells = swe_extra_params.h_ks_cells
            friction_factor_cells = swe_extra_params.friction_factor_cells
            Re_cells = swe_extra_params.Re_cells

            ManningN_cells_from_formula = deepcopy(ManningN_cells)

            #If the UDE choice is ManningN_h_Umag_ks, then compute the Reynolds number, h_ks, and friction factor
            if settings.UDE_settings.UDE_choice == "ManningN_h_Umag_ks"
                
                # Computer the Reynolds number

                ν = 1.0e-6 # kinematic viscosity of water (m^2/s)
                Re_cells = Umag .* h_array ./ ν

                Zygote.ignore() do
                    println("Re_cells: ", ForwardDiff.value.(Re_cells[1:5]))
                end

                # Compute h/ks 
                h_ks_cells = h_array ./ ks_cells

                #compute alpha
                alpha = 1.0 ./ (1.0 .+ (Re_cells ./ 850.0) .^ 9)

                #compute beta
                beta = 1.0 ./ (1.0 .+ (Re_cells ./ (h_ks_cells .* 160.0)) .^ 2)

                # Compute the friction factor f in parts 
                part1 = (Re_cells ./ 24.0) .^ alpha
                part2 = (1.8 .* log10.(Re_cells ./ 2.1)) .^ (2.0 .* (1.0 .- alpha) .* beta)
                part3 = (2.0 .* log10.(11.8 .* h_ks_cells)) .^ (2.0 .* (1.0 .- alpha) .* (1.0 .- beta))

                # Compute the friction factor f
                friction_factor_cells = 1.0 ./ (part1 .* part2 .* part3)

                # Compute Manning's n = sqrt(f/8.0) * h^(1/6) /sqrt(9.81)
                ManningN_cells_from_formula = sqrt.(friction_factor_cells ./ 8.0) .* h_array .^ (1.0 ./ 6.0) ./ sqrt.(9.81)
            end

            scalar_data = [h_array, q_x_array, q_y_array, zb_i, WSE_array, ManningN_cells, ManningN_cell_truth, ks_cells, h_ks_cells, friction_factor_cells, Re_cells, ManningN_cells_from_formula, WSE_truth, h_truth, u_truth, v_truth]
            scalar_names = ["h", "hu", "hv", "zb_cell", "WSE", "ManningN_cells", "ManningN_cell_truth", "ks", "h_ks", "friction_factor_from_formula", "Re", "ManningN_from_formula", "WSE_truth", "h_truth", "u_truth", "v_truth"]

            vtk_fileName = @sprintf("UDE_iteration_results_%04d.vtk", iter_number)                

            file_path = joinpath(case_path, vtk_fileName ) 
            export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, field_name, field_type, field_value, scalar_data, scalar_names, vector_data, vector_names)    
        end
    end

    #save everything (loss and ODE solution) to a JSON file
    println("       Saving UDE results to JSON file ...")
    save_file_name = "UDE_results_losses_and_ODE_solution.json"
    save_file_path = joinpath(case_path, save_file_name)
    open(save_file_path, "w") do io
        JSON3.pretty(io, Dict(
            "nIter" => nIter,
            "loss_totals" => loss_totals,
            "loss_preds" => loss_preds,
            "loss_pred_WSEs" => loss_pred_WSEs,
            "loss_pred_uvs" => loss_pred_uvs,
            "inverted_WSEs" => inverted_WSEs,
            "inverted_hs" => inverted_hs,
            "inverted_us" => inverted_us,
            "inverted_vs" => inverted_vs,
            "ManningN_cells" => ManningN_cells,
            "ManningN_cell_truth" => ManningN_cell_truth,
            "ks_cells" => swe_extra_params.ks_cells,
            "h_ks_cells" => swe_extra_params.h_ks_cells,
            "friction_factor_cells" => swe_extra_params.friction_factor_cells,
            "Re_cells" => swe_extra_params.Re_cells
        ))
    end
end

# Define the loss function
function compute_loss_UDE(ode_f, Q0, tspan, p, settings, my_mesh_2D, swe_2D_constants, observed_data, wstill, hstill, data_type)
    # Solve the ODE (forward pass)
    pred = UDE_forward_simulation(ode_f, Q0, tspan, p, settings, swe_2D_constants)

    #compute the loss
    # Ensure type stability in loss computation
    loss_total = zero(data_type)
    loss_pred = zero(data_type)
    loss_pred_WSE = zero(data_type)
    loss_pred_uv = zero(data_type)

    ϵ = sqrt(eps(data_type))

    if pred.retcode == SciMLBase.ReturnCode.Success
        WSE_truth = Vector{Float64}(observed_data["WSE_truth"])
        #h_truth = Vector{Float64}(observed_data["h_truth"])
        u_truth = Vector{Float64}(observed_data["u_truth"])
        v_truth = Vector{Float64}(observed_data["v_truth"])

        zb_cells_temp = observed_data["zb_cell_truth"]

        # Get min and max from truth data
        WSE_truth_min = minimum(WSE_truth)
        WSE_truth_max = maximum(WSE_truth)
        WSE_truth_range = WSE_truth_max - WSE_truth_min

        xi_pred = @view pred[1:my_mesh_2D.numOfCells, end]
        h_pred = xi_pred .+ hstill

        q_x_pred = @view pred[my_mesh_2D.numOfCells+1:2*my_mesh_2D.numOfCells, end]
        q_y_pred = @view pred[2*my_mesh_2D.numOfCells+1:3*my_mesh_2D.numOfCells, end]

        u_pred = @. q_x_pred / (h_pred + ϵ)
        v_pred = @. q_y_pred / (h_pred + ϵ)

        WSE_pred = h_pred .+ zb_cells_temp        

        # Normalize the WSE loss using truth data range
        loss_wse = (WSE_pred .- WSE_truth) ./ (WSE_truth_range .+ ϵ)       

        #loss for free surface elevation mismatch
        if settings.UDE_settings.UDE_bWSE_loss
            loss_pred_WSE = sum(abs2, loss_wse)/my_mesh_2D.numOfCells  #mean square error
        end

        #loss for velocity mismatch
        if settings.UDE_settings.UDE_b_uv_loss      #if also include u in the loss 
            # Compute velocity magnitude scale
            vel_mag_max = maximum(sqrt.(u_truth.^2 + v_truth.^2))
            vel_scale = max(vel_mag_max, eps(eltype(u_truth)))

            # Normalize using the characteristic velocity scale
            loss_u = (u_pred .- u_truth) ./ vel_scale
            loss_v = (v_pred .- v_truth) ./ vel_scale

            # Compute loss with optional weighting
            loss_pred_uv = (sum(abs2, loss_u) + sum(abs2, loss_v))/my_mesh_2D.numOfCells  #mean square error
        end

        #combined loss due to free surface elevation mismatch and velocity mismatch
        loss_pred = loss_pred_WSE + loss_pred_uv      #loss_pred is the same as loss_total (for now)
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

    return loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, pred
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
    elseif settings.UDE_settings.UDE_ode_solver == "Euler()"
        ode_solver = Euler()
    elseif settings.UDE_settings.UDE_ode_solver == "AB3()"
        ode_solver = AB3()
    elseif settings.UDE_settings.UDE_ode_solver == "RK4()"
        ode_solver = RK4()
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
