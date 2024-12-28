#This file is used to parse the control configuration file and create the settings for the forward simulation, inversion, and sensitivity analysis


# Define structs for different settings groups
struct TimeSettings
    bUse_srhhydro_time_settings::Bool
    tspan::Tuple{Float64,Float64}
    dt::Float64
end

struct ForwardSimulationSettings
    solver::String
    ode_solver::String
    ode_solver_adaptive::Bool
    ode_solver_b_jac_sparsity::Bool
    nSave::Int
    initial_condition_options::String
    initial_condition_file_name::String
    initial_condition_values_from_file::Union{Dict,Nothing}
    initial_condition_constant_values::Vector{Float64}
    save_file_name::String
    save_solution_truth_file_name::String
end

struct InversionSettings
    active_param_names::Vector{String}
    parameter_initial_values_options::String
    parameter_initial_values_file_name::String
    parameter_initial_values_from_file::Union{Dict,Nothing}
    zb_initial_values::Vector{Float64}
    ManningN_initial_values::Vector{Float64}
    inlet_discharge_initial_values::Vector{Float64}
    bInversion_slope_loss::Bool
    bInversion_u_loss::Bool
    lower_bound_zb::Float64
    upper_bound_zb::Float64
    optimizer::String
    learning_rate::Float64
    max_iterations::Int
    inversion_truth_file_name::String
    forward_simulation_initial_condition_options::String
    forward_simulation_initial_condition_file_name::String
    forward_simulation_initial_condition_constant_values::Vector{Float64}
    forward_simulation_initial_condition_values_from_file::Union{Dict,Nothing}
    ode_solver::String
    ode_solver_adaptive::Bool
    ode_solver_b_jac_sparsity::Bool
    ode_solver_nSave::Int
    save_file_name::String
    save_loss_history_file_name::String
end

struct ControlSettings
    bVerbose::Bool
    srhhydro_file_name::String
    bPerform_Forward_Simulation::Bool
    bPerform_Inversion::Bool
    bPerform_Sensitivity_Analysis::Bool
    time_settings::TimeSettings
    forward_settings::Union{ForwardSimulationSettings,Nothing}
    inversion_settings::Union{InversionSettings,Nothing}
end

# Function to parse control file and create settings
function parse_control_file(control_file::String)
    control_dict = JSON3.read(open(control_file), Dict)

    # Parse time settings
    time_settings = TimeSettings(
        control_dict["time_settings"]["bUse_srhhydro_time_settings"],
        Tuple(Float64.(control_dict["time_settings"]["tspan"])),
        Float64(control_dict["time_settings"]["dt"])
    )

    # Parse forward simulation settings if needed
    forward_settings = nothing
    if control_dict["control_variables"]["bPerform_Forward_Simulation"]
        forward_settings = parse_forward_settings(control_dict["forward_simulation_options"])
    end

    # Parse inversion settings if needed
    inversion_settings = nothing
    if control_dict["control_variables"]["bPerform_Inversion"]
        inversion_settings = parse_inversion_settings(control_dict["inversion_options"])
    end

    # Create main control settings
    settings = ControlSettings(
        control_dict["bVerbose"],
        control_dict["control_variables"]["srhhydro_file_name"],
        control_dict["control_variables"]["bPerform_Forward_Simulation"],
        control_dict["control_variables"]["bPerform_Inversion"],
        control_dict["control_variables"]["bPerform_Sensitivity_Analysis"],
        time_settings,
        forward_settings,
        inversion_settings
    )

    # Validate settings
    if settings.bPerform_Forward_Simulation + settings.bPerform_Inversion + settings.bPerform_Sensitivity_Analysis != 1
        error("Currently only one at a time: forward simulation, inversion, or sensitivity analysis.")
    end

    # Use settings throughout the code
    if settings.bVerbose
        println("--------------------------------")
        println("Control variables:")
        println("   srhhydro_file_name = ", settings.srhhydro_file_name)
        println("   bPerform_Forward_Simulation = ", settings.bPerform_Forward_Simulation)
        println("   bPerform_Inversion = ", settings.bPerform_Inversion)
        println("   bPerform_Sensitivity_Analysis = ", settings.bPerform_Sensitivity_Analysis)

        println("Time settings (overwrite SRH-2D case if bUse_srhhydro_time_settings is false):")
        println("    bUse_srhhydro_time_settings = ", settings.time_settings.bUse_srhhydro_time_settings)
        println("    tspan = ", settings.time_settings.tspan)
        println("    dt = ", settings.time_settings.dt)

        if settings.bPerform_Forward_Simulation
            println("Forward simulation is to be performed.")
            println("Forward simulation options:")
            println("    solver = ", settings.forward_settings.solver)
            println("    ode_solver = ", settings.forward_settings.ode_solver)
            println("    ode_solver_adaptive = ", settings.forward_settings.ode_solver_adaptive)
            println("    ode_solver_b_jac_sparsity = ", settings.forward_settings.ode_solver_b_jac_sparsity)
            println("    nSave = ", settings.forward_settings.nSave)
            println("    initial_condition_options = ", settings.forward_settings.initial_condition_options)
            println("    initial_condition_file_name = ", settings.forward_settings.initial_condition_file_name)
            println("    save_file_name = ", settings.forward_settings.save_file_name)
            println("    save_solution_truth_file_name = ", settings.forward_settings.save_solution_truth_file_name)
        else
            println("No forward simulation is to be performed.")
        end

        if settings.bPerform_Inversion
            println("Inversion is to be performed.")
            println("Inversion options:")
            println("    active_param_names = ", settings.inversion_settings.active_param_names)
            println("    inversion_parameter_initial_values_options = ", settings.inversion_settings.parameter_initial_values_options)
            println("    inversion_parameter_initial_values_file_name = ", settings.inversion_settings.parameter_initial_values_file_name)

            if settings.inversion_settings.parameter_initial_values_options == "from_file"
                println("    inversion_parameter_initial_values_from_file = ", settings.inversion_settings.parameter_initial_values_from_file)
            elseif settings.inversion_settings.parameter_initial_values_options == "constant"
                println("    inversion_zb_initial_values = ", settings.inversion_settings.zb_initial_values)
                println("    inversion_ManningN_initial_values = ", settings.inversion_settings.ManningN_initial_values)
                println("    inversion_inlet_discharge_initial_values = ", settings.inversion_settings.inlet_discharge_initial_values)
            else
                error("Invalid inversion_parameter_initial_values_options: $inversion_parameter_initial_values_options. Supported options: from_file, constant.")
            end

            println("    inversion_bInversion_slope_loss = ", settings.inversion_settings.bInversion_slope_loss)
            println("    inversion_bInversion_u_loss = ", settings.inversion_settings.bInversion_u_loss)
            println("    inversion_lower_bound_zb = ", settings.inversion_settings.lower_bound_zb)
            println("    inversion_upper_bound_zb = ", settings.inversion_settings.upper_bound_zb)
            println("    optimizer = ", settings.inversion_settings.optimizer)
            println("    learning_rate = ", settings.inversion_settings.learning_rate)
            println("    max_iterations = ", settings.inversion_settings.max_iterations)
            println("    inversion_truth_file_name = ", settings.inversion_settings.inversion_truth_file_name)
            println("    inversion_forward_simulation_initial_condition_options = ", settings.inversion_settings.forward_simulation_initial_condition_options)
            println("    inversion_forward_simulation_initial_condition_file_name = ", settings.inversion_settings.forward_simulation_initial_condition_file_name)
            println("    inversion_forward_simulation_initial_condition_constant_values = ", settings.inversion_settings.forward_simulation_initial_condition_constant_values)
            println("    ode_solver = ", settings.inversion_settings.ode_solver)
            println("    ode_solver_adaptive = ", settings.inversion_settings.ode_solver_adaptive)
            println("    ode_solver_b_jac_sparsity = ", settings.inversion_settings.ode_solver_b_jac_sparsity)
        else
            println("No inversion is to be performed.")
        end

        if settings.bPerform_Sensitivity_Analysis
            println("Sensitivity analysis is to be performed.")
        else
            println("No sensitivity analysis is to be performed.")
        end

        println("--------------------------------")
    end

    return settings
end

# Helper functions to parse specific settings
function parse_forward_settings(options::Dict)
    # Handle initial condition values from file if specified
    initial_condition_values = nothing
    if options["forward_simulation_initial_condition_options"] == "from_file"
        file_path = joinpath(dirname(@__FILE__),
            options["forward_simulation_initial_condition_file_name"])
        initial_condition_values = JSON3.read(open(file_path), Dict)
    end

    ForwardSimulationSettings(
        options["forward_simulation_solver"],                    # solver
        options["forward_simulation_ode_solver"],               # ode_solver
        options["forward_simulation_adaptive"],                 # adaptive
        options["forward_simulation_b_jac_sparsity"],          # b_jac_sparsity
        Int(options["forward_simulation_nSave"]),              # nSave
        options["forward_simulation_initial_condition_options"],  # initial_condition_options
        options["forward_simulation_initial_condition_file_name"], # initial_condition_file_name
        initial_condition_values,                               # initial_condition_values_from_file
        Float64.(options["forward_simulation_initial_condition_constant_values"]), # initial_condition_constant_values
        options["forward_simulation_save_file_name"],           # save_file_name
        options["forward_simulation_save_solution_truth_file_name"] # save_solution_truth_file_name
    )
end

function parse_inversion_settings(options::Dict)
    # Handle parameter initial values from file if specified
    parameter_values_from_file = nothing
    if options["inversion_parameter_initial_values_options"] == "from_file"
        file_path = joinpath(dirname(@__FILE__),
            options["inversion_parameter_initial_values_file_name"])
        parameter_values_from_file = JSON3.read(open(file_path), Dict)
    end

    forward_simulation_initial_condition_values_from_file = nothing
    if options["inversion_forward_simulation_initial_condition_options"] == "from_file"
        file_path = joinpath(dirname(@__FILE__),
            options["inversion_forward_simulation_initial_condition_file_name"])
        forward_simulation_initial_condition_values_from_file = JSON3.read(open(file_path), Dict)
    end

    InversionSettings(
        String.(options["active_param_names"]),                # active_param_names
        options["inversion_parameter_initial_values_options"],  # parameter_initial_values_options
        options["inversion_parameter_initial_values_file_name"], # parameter_initial_values_file_name
        parameter_values_from_file,                            # parameter_initial_values_from_file
        Float64.(options["inversion_zb_initial_values"]),      # zb_initial_values
        Float64.(options["inversion_ManningN_initial_values"]), # ManningN_initial_values
        Float64.(options["inversion_inlet_discharge_initial_values"]), # inlet_discharge_initial_values
        options["inversion_bInversion_slope_loss"],            # bInversion_slope_loss
        options["inversion_bInversion_u_loss"],                # bInversion_u_loss
        Float64(options["inversion_lower_bound_zb"]),          # lower_bound_zb
        Float64(options["inversion_upper_bound_zb"]),          # upper_bound_zb
        options["inversion_optimizer"],                        # optimizer
        Float64(options["inversion_learning_rate"]),           # learning_rate
        Int(options["inversion_max_iterations"]),              # max_iterations
        options["inversion_truth_file_name"],         # inversion_truth_file_name
        options["inversion_forward_simulation_initial_condition_options"], # forward_simulation_initial_condition_options
        options["inversion_forward_simulation_initial_condition_file_name"], # forward_simulation_initial_condition_file_name
        Float64.(options["inversion_forward_simulation_initial_condition_constant_values"]), # forward_simulation_initial_condition_constant_values,
        forward_simulation_initial_condition_values_from_file, # forward_simulation_initial_condition_values_from_file
        options["inversion_ode_solver_options"]["ode_solver"], # ode_solver
        options["inversion_ode_solver_options"]["ode_solver_adaptive"], # ode_solver_adaptive
        options["inversion_ode_solver_options"]["ode_solver_b_jac_sparsity"], # ode_solver_b_jac_sparsity
        Int(options["inversion_ode_solver_options"]["ode_solver_nSave"]), # ode_solver_nSave
        options["inversion_save_file_name"],                   # save_file_name
        get(options, "inversion_save_loss_history_file_name", "loss_history.jld2") # save_loss_history_file_name
    )
end