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

    # Create ODEProblem (params_vector is the initial guess for the parameters; it will be updated by the optimizer)       
    #ode_prob = ODEProblem(ode_f, Q0, swe_2D_constants.tspan, params_vector)

    #perform the inversion
    println("   Performing inversion ...\n")
    #sol, ITER, LOSS, PARS = optimize_parameters_inversion(ode_f, Q0, params_vector, settings, my_mesh_2D, swe_2D_constants, observed_data, active_param_name, case_path)

    #compute all the losses and ODE solutions for each inversion iteration
    println("   Computing all the losses and ODE solutions for each inversion iteration ...")
    compute_all_losses_inversion(ode_f, Q0, settings, nodeCoordinates, my_mesh_2D, swe_2D_constants, observed_data, active_param_name, case_path)

    #save the inversion results
    #jldsave(joinpath(case_path, settings.inversion_settings.save_file_name); ITER, LOSS, PARS)

    #process the inversion results
    println("   Post-processing inversion results ...")

    #process inversion results
    #Hydrograd.postprocess_inversion_results_swe_2D(settings, my_mesh_2D, nodeCoordinates, zb_cell_truth, h_truth, u_truth, v_truth, WSE_truth, case_path)

end

# compute all the losses for each inversion iteration
function compute_all_losses_inversion(ode_f::ODEFunction, Q0::AbstractVector{T1}, settings::ControlSettings, nodeCoordinates::Matrix{T2}, 
    my_mesh_2D::mesh_2D, swe_2D_constants::swe_2D_consts, observed_data::Dict{String,Vector{T2}}, active_param_name::String, case_path::String) where {T1<:Real,T2<:Real}

    data_type = eltype(Q0)
    
    #find the list of files of "inversion_callback_save_iter_<iter_number>.json"
    iter_files = filter(file -> occursin("inversion_callback_save_iter_", file) && endswith(file, ".json"), readdir(case_path))

    #println("iter_files = ", iter_files)

    #sort the files by the iteration number
    iter_files = sort(iter_files, by = file -> begin
        # Extract the number between "iter_" and ".json"
        m = match(r"inversion_callback_save_iter_(\d+)\.json$", file)
        if isnothing(m)
            error("Invalid filename format: $file")
        end
        parse(Int, m.captures[1])
    end)

    #print the total number of iterations
    nIter = length(iter_files)
    println("       Total number of inversion iterations: ", nIter)

    #get the truth data
    WSE_truth = observed_data["WSE_truth"]
    h_truth = observed_data["h_truth"]
    u_truth = observed_data["u_truth"]
    v_truth = observed_data["v_truth"]
    zb_cell_truth = observed_data["zb_cell_truth"]
    ManningN_cell_truth = observed_data["ManningN_cell_truth"]
    inlet_discharges_truth = observed_data["inlet_discharges_truth"]    


    #initialize the losses
    loss_totals = []
    loss_preds = []
    loss_pred_WSEs = []
    loss_pred_uvs = []
    loss_bounds = []
    loss_slopes = []
    inverted_WSEs = []
    inverted_hs = []
    inverted_us = []
    inverted_vs = []
    
    #load the files
    for file in iter_files
        println("       Loading inversion results from file: ", file)
        data = JSON3.read(open(joinpath(case_path, file)), Dict{String,Any})

        iter_number = parse(Int, match(r"inversion_callback_save_iter_(\d+)\.json$", file).captures[1])
        iter_number_from_file = data["iter_number"]
        @assert iter_number == iter_number_from_file "The iteration number in the file name and the file content do not match: iter_number = $iter_number, iter_number_from_file = $iter_number_from_file"

        θ_from_file = convert(Vector{Float64}, data["theta"])
        total_loss_from_file = convert(Float64, data["loss_total"])

        #@show typeof(θ_from_file)

        loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_bound, loss_slope, ode_pred = compute_loss_inversion(ode_f, Q0, θ_from_file, settings,
            my_mesh_2D, swe_2D_constants, observed_data, active_param_name, data_type)
            

        #assert the total loss is the same as the one in the file
        @assert abs(loss_total - total_loss_from_file) < 1e-3 "The total loss in the file and the computed loss do not match: loss_total = $loss_total, total_loss_from_file = $total_loss_from_file"

        #append the losses to the accumulators
        push!(loss_totals, loss_total)
        push!(loss_preds, loss_pred)
        push!(loss_pred_WSEs, loss_pred_WSE)
        push!(loss_pred_uvs, loss_pred_uv)
        push!(loss_bounds, loss_bound)
        push!(loss_slopes, loss_slope)

        #if the inversion is for zb, then curPars is the zb_cells_param. Otherwise, curPar is something else such as the Manning's n or the inlet discharges. 
        if settings.inversion_settings.active_param_names == ["zb"]
            zb_i = θ_from_file

            #compute S0
            #_, _, S0_i = interploate_zb_from_cell_to_face_and_compute_S0(my_mesh_2D, zb_i)
        else
            zb_i = zb_cell_truth
        end

        @show size(ode_pred.u[end])
        @show size(zb_i)

        #append the ODE solution to the accumulator
        push!(inverted_WSEs, ode_pred.u[end][1:my_mesh_2D.numOfCells] .+ zb_i)
        push!(inverted_hs, ode_pred.u[end][1:my_mesh_2D.numOfCells])
        push!(inverted_us, ode_pred.u[end][my_mesh_2D.numOfCells+1:2*my_mesh_2D.numOfCells])
        push!(inverted_vs, ode_pred.u[end][2*my_mesh_2D.numOfCells+1:3*my_mesh_2D.numOfCells])


        #save the inverted results to vtk
        bSave_vtk = true
        if bSave_vtk
            println("       Saving forward simulation results to vtk file for iteration: ", iter_number)

            field_name = "iter_number"
            field_type = "integer"
            field_value = iter_number

            h_array = ode_pred.u[end][1:my_mesh_2D.numOfCells]
            q_x_array = ode_pred.u[end][my_mesh_2D.numOfCells+1:2*my_mesh_2D.numOfCells]
            q_y_array = ode_pred.u[end][2*my_mesh_2D.numOfCells+1:3*my_mesh_2D.numOfCells]

            u_temp = q_x_array ./ h_array
            v_temp = q_y_array ./ h_array
            U_vector = hcat(u_temp, v_temp)
                            
            vector_data = [U_vector] 
            vector_names = ["U"]

            WSE_array = h_array + zb_i
                
            scalar_data = [h_array, q_x_array, q_y_array, zb_i, WSE_array]
            scalar_names = ["h", "hu", "hv", "zb_cell", "WSE"]

            vtk_fileName = @sprintf("forward_simulation_results_%04d.vtk", iter_number)
                
            file_path = joinpath(case_path, vtk_fileName ) 
            export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, field_name, field_type, field_value, scalar_data, scalar_names, vector_data, vector_names)    
        end
    end

    #save everything (loss and ODE solution) to a JSON file
    println("       Saving inversion results to JSON file ...")
    save_file_name = "inversion_results_losses_and_ODE_solution.json"
    save_file_path = joinpath(case_path, save_file_name)
    open(save_file_path, "w") do io
        JSON3.pretty(io, Dict(
            "nIter" => nIter,
            "loss_totals" => loss_totals,
            "loss_preds" => loss_preds,
            "loss_pred_WSEs" => loss_pred_WSEs,
            "loss_pred_uvs" => loss_pred_uvs,
            "loss_bounds" => loss_bounds,
            "loss_slopes" => loss_slopes,
            "inverted_WSEs" => inverted_WSEs,
            "inverted_hs" => inverted_hs,
            "inverted_us" => inverted_us,
            "inverted_vs" => inverted_vs
        ))
    end
    

end

# Define the loss function
function compute_loss_inversion(ode_f::ODEFunction, Q0::AbstractVector{T1}, p::AbstractVector{T2}, settings::ControlSettings, 
    my_mesh_2D::mesh_2D, swe_2D_constants::swe_2D_consts, observed_data::Dict{String,Vector{T3}}, active_param_name::String, data_type::DataType) where {T1<:Real,T2<:Real,T3<:Real}

    # Solve the ODE (forward pass)

    # Create ODEProblem        
    ode_prob = ODEProblem(ode_f, Q0, swe_2D_constants.tspan, p)

    # Create ODE solver
    if settings.inversion_settings.ode_solver == "Tsit5()"
        ode_solver = Tsit5()
    elseif settings.inversion_settings.ode_solver == "TRBDF2()"
        ode_solver = TRBDF2()
    elseif settings.inversion_settings.ode_solver == "Rosenbrock23()"
        ode_solver = Rosenbrock23()
    elseif settings.inversion_settings.ode_solver == "Euler()"
        ode_solver = Euler()
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
    if settings.inversion_settings.ode_solver_sensealg == "ZygoteAdjoint()"  #fails
        ode_solver_sensealg = ZygoteAdjoint()
    elseif settings.inversion_settings.ode_solver_sensealg == "InterpolatingAdjoint(autojacvec=nothing)"
        ode_solver_sensealg = InterpolatingAdjoint(autojacvec=nothing)
    elseif settings.inversion_settings.ode_solver_sensealg == "ForwardDiffSensitivity()"
        ode_solver_sensealg = ForwardDiffSensitivity()
    elseif settings.inversion_settings.ode_solver_sensealg == "BacksolveAdjoint(autojacvec=nothing)"
        ode_solver_sensealg = BacksolveAdjoint(autojacvec=nothing)    
    elseif settings.inversion_settings.ode_solver_sensealg == "GaussAdjoint(autojacvec=ZygoteVJP())"
        ode_solver_sensealg = GaussAdjoint(autojacvec=ZygoteVJP())
    elseif settings.inversion_settings.ode_solver_sensealg == "QuadratureAdjoint(autojacvec=nothing)"
        ode_solver_sensealg = QuadratureAdjoint(autojacvec = nothing)
    elseif settings.inversion_settings.ode_solver_sensealg == "ReverseDiffAdjoint()"
        ode_solver_sensealg = ReverseDiffAdjoint()
    elseif settings.inversion_settings.ode_solver_sensealg == "ForwardSensitivity()"
        ode_solver_sensealg = ForwardSensitivity()
    elseif settings.inversion_settings.ode_solver_sensealg == "AutoEnzyme()"
        ode_solver_sensealg = AutoEnzyme()                
    else
        error("Not implemented yet")
    end

    #@show ode_solver_sensealg

    #define the time for saving the results for the ODE solver
    dt_save = (swe_2D_constants.tspan[2] - swe_2D_constants.tspan[1]) / settings.inversion_settings.ode_solver_nSave
    t_save = swe_2D_constants.tspan[1]:dt_save:swe_2D_constants.tspan[2]

    # Solve the ODE with the updated parameters p
    #@time 
    ode_pred = solve(ode_prob, ode_solver, p=p, adaptive=settings.inversion_settings.ode_solver_adaptive, dt=swe_2D_constants.dt, saveat=t_save, sensealg=ode_solver_sensealg)

    Zygote.ignore() do
        if settings.bVerbose
            #@show ode_pred.retcode
            #@show typeof(ode_pred)
            #@show size(ode_pred)
            #@show ode_pred
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

    if ode_pred.retcode == SciMLBase.ReturnCode.Success
        WSE_truth = Vector{Float64}(observed_data["WSE_truth"])
        #h_truth = Vector{Float64}(observed_data["h_truth"])
        u_truth = Vector{Float64}(observed_data["u_truth"])
        v_truth = Vector{Float64}(observed_data["v_truth"])

        #Get the bed elevation at cells: if active_param_name == "zb", then use the initial values; otherwise, use the truth data
        if active_param_name == "zb"
            zb_cells_temp = p
        else
            zb_cells_temp = observed_data["zb_cell_truth"]
        end

        # Add small epsilon to prevent division by zero
        ϵ = sqrt(eps(data_type))

        Zygote.ignore() do
            #@show size(ode_pred)
        end

        h_pred = @view ode_pred[1:my_mesh_2D.numOfCells, end]
        q_x_pred = @view ode_pred[my_mesh_2D.numOfCells+1:2*my_mesh_2D.numOfCells, end]
        q_y_pred = @view ode_pred[2*my_mesh_2D.numOfCells+1:3*my_mesh_2D.numOfCells, end]
        u_pred = @. q_x_pred / (h_pred + ϵ)
        v_pred = @. q_y_pred / (h_pred + ϵ)

        WSE_pred = h_pred .+ zb_cells_temp

        # Get min and max from truth data
        WSE_truth_min = minimum(WSE_truth)
        WSE_truth_max = maximum(WSE_truth)
        WSE_truth_range = WSE_truth_max - WSE_truth_min

        # Normalize the WSE loss using truth data range
        loss_wse = (WSE_pred .- WSE_truth) ./ (WSE_truth_range .+ ϵ)

        #loss for free surface elevation mismatch
        if settings.inversion_settings.bInversion_WSE_loss
            loss_pred_WSE = sum(abs2, loss_wse)/my_mesh_2D.numOfCells  #mean square error
        end

        #loss for velocity mismatch
        if settings.inversion_settings.bInversion_uv_loss      #if also include velocity mismatch in the loss 
            # Get characteristic velocity scale from truth data
            #u_max = maximum(abs.(u_truth))
            #v_max = maximum(abs.(v_truth))
            #vel_scale = max(u_max, v_max, eps(eltype(u_truth)))  # Use larger component

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
        loss_pred = loss_pred_WSE + loss_pred_uv

        #loss for parameter bound regularization
        if settings.inversion_settings.bInversion_bound_loss
            loss_bound = compute_bound_loss(p, settings.inversion_settings.parameter_value_lower_bound, settings.inversion_settings.parameter_value_upper_bound)/my_mesh_2D.numOfCells  #mean square error
        end

        #loss for bed slope regularization
        if settings.inversion_settings.bInversion_slope_loss    #if bed slope is included in the loss 
            loss_slope = calc_slope_loss(zb_cells_temp, settings, my_mesh_2D)  #slope loss is already normalized by the number of cells
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

    return loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_bound, loss_slope, ode_pred
end

function optimize_parameters_inversion(ode_f, Q0, p_init, settings, my_mesh_2D, swe_2D_constants, observed_data, active_param_name, case_path)
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

        loss_total, loss_pred, loss_pred_WSE, loss_pred_uv, loss_bound, loss_slope, ode_pred = compute_loss_inversion(ode_f, Q0, θ, settings,
            my_mesh_2D, swe_2D_constants, observed_data, active_param_name, data_type)

        return loss_total
    end

    # Define AD type choice for optimization's gradient computation: https://docs.sciml.ai/Optimization/stable/API/ad/
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
    elseif settings.inversion_settings.inversion_sensealg == "AutoEnzyme()"
        adtype = Optimization.AutoEnzyme()
    elseif settings.inversion_settings.inversion_sensealg == "AutoReverseDiff()"
        adtype = Optimization.AutoReverseDiff(compile=false)
    elseif settings.inversion_settings.inversion_sensealg == "AutoForwardDiff()"
        adtype = Optimization.AutoForwardDiff()
    elseif settings.inversion_settings.inversion_sensealg == "ForwardSensitivity()"
        adtype = Optimization.ForwardSensitivity()
    elseif settings.inversion_settings.inversion_sensealg == "AutoFiniteDiff()"
        adtype = Optimization.AutoFiniteDiff()
    elseif settings.inversion_settings.inversion_sensealg == "AutoModelingToolkit()"
        adtype = Optimization.AutoModelingToolkit()
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

    # Define the bounds for the parameter (only applicable for some optimizers which support lb and ub such as LBFGS)
    lb_p = ones(my_mesh_2D.numOfCells)
    ub_p = ones(my_mesh_2D.numOfCells)

    lb_p .= settings.inversion_settings.parameter_value_lower_bound
    ub_p .= settings.inversion_settings.parameter_value_upper_bound

    # Define the bounds for the parameter (only applicable for some optimizers which support lb and ub such as LBFGS)
    # Note: lb and ub are not used for the optimizer currently
    #lb_p = zeros(my_mesh_2D.numOfCells)
    #lb_p .= -0.1
    #ub_p = zeros(my_mesh_2D.numOfCells)
    #ub_p .= 0.3

    #From SciMLSensitivity documentation: https://docs.sciml.ai/Optimization/stable/API/optimization_problem/
    #OptimizationProblem{iip}(f, u0, p = SciMLBase.NullParameters(),;
    #                          lb = nothing,
    #                          ub = nothing,
    #                          lcons = nothing,
    #                          ucons = nothing,
    #                          sense = nothing,
    #                          kwargs...)
    # optprob = nothing
    # if settings.inversion_settings.inversion_optimizers[1] == "Descent"  #Descent does not support lb and ub
    #     optprob = OptimizationProblem(optf, p_init) 
    # elseif settings.inversion_settings.inversion_optimizers[1] == "Adam"
    #     optprob = OptimizationProblem(optf, p_init, lb=lb_p, ub=ub_p)  
    # elseif settings.inversion_settings.inversion_optimizers[1] == "LBFGS" #LBFGS does support lb and ub
    #     optprob = OptimizationProblem(optf, p_init, lb=lb_p, ub=ub_p)  
    # else
    #     error("Not implemented yet")
    # end

    optprob = OptimizationProblem(optf, p_init) 

    #define the accumulators for the inversion results
    ITER = []        #iteration number accumulator
    LOSS = []        # Loss accumulator
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

        #append the loss and parameters to the accumulators
        append!(ITER, iter_number)
        append!(LOSS, [loss_total])
        append!(PARS, [θ])

        Zygote.ignore() do
            println("       iter_number = ", iter_number, ", loss_total = ", loss_total, ", inversion_elapsed_seconds = ", inversion_elapsed_seconds)
            #println("       optimizer_state = ", optimizer_state)

            #if not the first iteration, compute and print the max and min of the parameter change from the previous iteration
            if iter_number > 1
                param_change = θ - PARS[end-1]
                println("           max(param_change) = ", maximum(param_change), ", min(param_change) = ", minimum(param_change))
            end

            #save the inversion results (for just the current iteration; insurance in case the inversion is interrupted)
            if iter_number % settings.inversion_settings.save_frequency == 0
                #jldsave(joinpath(case_path, "inversion_callback_save_iter_$(iter_number).jld2"); iter_number, loss_total, θ)

                call_back_save_dict = Dict("iter_number" => iter_number, "loss_total" => loss_total, "theta" => θ)
                open(joinpath(case_path, "inversion_callback_save_iter_$(iter_number).json"), "w") do io
                    JSON3.pretty(io, call_back_save_dict)
                end
            end

            #checkpoint the inversion results (in case the inversion is interrupted)
            if settings.inversion_settings.save_checkpoint && iter_number % settings.inversion_settings.checkpoint_frequency == 0
                jldsave(joinpath(case_path, "checkpoint_inversion_iter_$(iter_number).jld2"); ITER, LOSS, PARS)
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
        if optimizer_choice == "Descent"
            optimizer = OptimizationOptimisers.Descent(settings.inversion_settings.inversion_learning_rates[iOptimizer])
        elseif optimizer_choice == "Adam"
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

        if optimizer_choice == "Descent" 
            sol = solve(optprob, optimizer, callback=callback, maxiters=max_iterations,
                abstol=settings.inversion_settings.inversion_abs_tols[iOptimizer], reltol=settings.inversion_settings.inversion_rel_tols[iOptimizer])
        elseif optimizer_choice == "Adam" 
            #@time 
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

    return sol, ITER, LOSS, PARS
end

# Bound loss function compatible with both Zygote and ForwardDiff
function compute_bound_loss(params::Vector{T}, lower_bound::Real, upper_bound::Real) where {T}
    loss = zero(T)
    bound_range = upper_bound - lower_bound
    
    for p in params
        lower_violation = max(zero(T), lower_bound - p) / bound_range
        upper_violation = max(zero(T), p - upper_bound) / bound_range
        loss += sum(abs2, lower_violation) + sum(abs2, upper_violation)
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
    S0 = compute_scalar_gradients_from_cells(my_mesh_2D, params)

    #compute the magnitude of S0
    S0_mag = sqrt.(S0[:, 1] .^ 2 + S0[:, 2] .^ 2 .+ eps(T))

    #compute the loss due to slope regularization
    loss = compute_bound_loss(S0_mag, 0.0, settings.inversion_settings.parameter_slope_upper_bound)

    #normalize the loss by the number of cells
    loss = loss / my_mesh_2D.numOfCells

    return loss
end