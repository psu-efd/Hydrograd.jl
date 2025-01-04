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
    bInversion_WSE_loss::Bool    
    bInversion_uv_loss::Bool    
    bInversion_bound_loss::Bool
    bInversion_slope_loss::Bool
    parameter_value_lower_bound::Float64
    parameter_value_upper_bound::Float64
    parameter_slope_upper_bound::Float64
    optimizer::String
    learning_rate::Float64
    max_iterations::Int
    save_frequency::Int
    save_checkpoint::Bool
    checkpoint_frequency::Int
    inversion_sensealg::String
    inversion_truth_file_name::String
    forward_simulation_initial_condition_options::String
    forward_simulation_initial_condition_file_name::String
    forward_simulation_initial_condition_constant_values::Vector{Float64}
    forward_simulation_initial_condition_values_from_file::Union{Dict,Nothing}
    ode_solver::String
    ode_solver_adaptive::Bool
    ode_solver_sensealg::String
    ode_solver_b_jac_sparsity::Bool
    ode_solver_nSave::Int
    save_file_name::String
    save_loss_history_file_name::String
    save_parameters_history_file_name::String
end

struct SensitivityAnalysisSettings
    active_param_names::Vector{String}
    parameter_values_options::String
    parameter_values_file_name::String
    parameter_values_from_file::Union{Dict,Nothing}
    zb_values::Vector{Float64}
    ManningN_values::Vector{Float64}
    inlet_discharge_values::Vector{Float64}
    forward_simulation_initial_condition_options::String
    forward_simulation_initial_condition_file_name::String
    forward_simulation_initial_condition_constant_values::Vector{Float64}
    forward_simulation_initial_condition_values_from_file::Union{Dict,Nothing}
    ode_solver::String
    ode_solver_adaptive::Bool
    ode_solver_b_jac_sparsity::Bool
    ode_solver_nSave::Int
    save_file_name::String
end

struct UDESettings
    UDE_mode::String
    UDE_choice::String
    UDE_NN_config::Dict
    UDE_bWSE_loss::Bool
    UDE_b_uv_loss::Bool
    UDE_optimizer::String
    UDE_learning_rate::Float64
    UDE_max_iterations::Int
    UDE_training_save_frequency::Int
    UDE_training_save_checkpoint::Bool
    UDE_training_checkpoint_frequency::Int
    UDE_sensealg::String
    UDE_truth_file_name::String
    forward_simulation_initial_condition_options::String
    forward_simulation_initial_condition_file_name::String
    forward_simulation_initial_condition_constant_values::Vector{Float64}
    forward_simulation_initial_condition_values_from_file::Union{Dict,Nothing}
    UDE_ode_solver::String
    UDE_ode_solver_adaptive::Bool
    UDE_ode_solver_sensealg::String
    UDE_ode_solver_b_jac_sparsity::Bool
    UDE_ode_solver_nSave::Int
    UDE_training_save_file_name::String
    UDE_training_save_loss_history_file_name::String
    UDE_training_save_parameters_history_file_name::String
    UDE_training_save_NN_weights_state_file_name::String
    UDE_inference_save_file_name::String
end


struct ControlSettings
    bVerbose::Bool
    srhhydro_file_name::String
    bPerform_Forward_Simulation::Bool
    bPerform_Inversion::Bool
    bPerform_Sensitivity_Analysis::Bool
    bPerform_UDE::Bool
    time_settings::TimeSettings
    forward_settings::Union{ForwardSimulationSettings,Nothing}
    inversion_settings::Union{InversionSettings,Nothing}
    sensitivity_analysis_settings::Union{SensitivityAnalysisSettings,Nothing}
    UDE_settings::Union{UDESettings,Nothing}
end

# Function to parse control file and create settings
function parse_control_file(control_file::String)
    control_dict = JSON3.read(open(control_file), Dict)
    control_file_dir = dirname(control_file)  # Get directory of control file

    # check whether any of the task flags are present
    if !(haskey(control_dict["control_variables"], "bPerform_Forward_Simulation") &&
         haskey(control_dict["control_variables"], "bPerform_Inversion") &&
         haskey(control_dict["control_variables"], "bPerform_Sensitivity_Analysis") &&
         haskey(control_dict["control_variables"], "bPerform_UDE"))
        error("None of the task flags (bPerform_Forward_Simulation, bPerform_Inversion, bPerform_Sensitivity_Analysis, and bPerform_UDE) are present. Please set at least one of the task flags to true.")
    end

    if !any([control_dict["control_variables"]["bPerform_Forward_Simulation"],
            control_dict["control_variables"]["bPerform_Inversion"],
            control_dict["control_variables"]["bPerform_Sensitivity_Analysis"],
            control_dict["control_variables"]["bPerform_UDE"]])
        error("No task flag is set true. Please set at least one of the task flags to true.")
    end

    # Parse time settings (always required)
    time_settings = TimeSettings(
        control_dict["time_settings"]["bUse_srhhydro_time_settings"],
        Tuple(Float64.(control_dict["time_settings"]["tspan"])),
        Float64(control_dict["time_settings"]["dt"])
    )

    # Parse forward simulation settings if needed
    forward_settings = nothing
    if control_dict["control_variables"]["bPerform_Forward_Simulation"]
        forward_settings = parse_forward_settings(control_dict["forward_simulation_options"], control_file_dir)
    end

    # Parse inversion settings if needed
    inversion_settings = nothing
    if control_dict["control_variables"]["bPerform_Inversion"]
        inversion_settings = parse_inversion_settings(control_dict["inversion_options"], control_file_dir)
    end

    # Parse sensitivity analysis settings if needed
    sensitivity_analysis_settings = nothing
    if control_dict["control_variables"]["bPerform_Sensitivity_Analysis"]
        sensitivity_analysis_settings = parse_sensitivity_analysis_settings(control_dict["sensitivity_analysis_options"], control_file_dir)
    end

    # Parse UDE settings if needed
    UDE_settings = nothing
    if control_dict["control_variables"]["bPerform_UDE"]
        UDE_settings = parse_UDE_settings(control_dict["UDE_options"], control_file_dir)
    end

    # Create main control settings
    settings = ControlSettings(
        control_dict["bVerbose"],
        control_dict["control_variables"]["srhhydro_file_name"],
        control_dict["control_variables"]["bPerform_Forward_Simulation"],
        control_dict["control_variables"]["bPerform_Inversion"],
        control_dict["control_variables"]["bPerform_Sensitivity_Analysis"],
        control_dict["control_variables"]["bPerform_UDE"],
        time_settings,
        forward_settings,
        inversion_settings,
        sensitivity_analysis_settings,
        UDE_settings
    )

    # Validate settings
    if settings.bPerform_Forward_Simulation + settings.bPerform_Inversion + settings.bPerform_Sensitivity_Analysis + settings.bPerform_UDE != 1
        error("Currently only one at a time: forward simulation, inversion, sensitivity analysis, or UDE.")
    end

    # Use settings throughout the code
    if settings.bVerbose
        println("--------------------------------")
        println("Control variables:")
        println("   srhhydro_file_name = ", settings.srhhydro_file_name)
        println("   bPerform_Forward_Simulation = ", settings.bPerform_Forward_Simulation)
        println("   bPerform_Inversion = ", settings.bPerform_Inversion)
        println("   bPerform_Sensitivity_Analysis = ", settings.bPerform_Sensitivity_Analysis)
        println("   bPerform_UDE = ", settings.bPerform_UDE)

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
                error("Invalid inversion_parameter_initial_values_options: $(settings.inversion_settings.parameter_initial_values_options). Supported options: from_file, constant.")
            end

            println("    inversion_bInversion_WSE_loss = ", settings.inversion_settings.bInversion_WSE_loss)
            println("    inversion_bInversion_uv_loss = ", settings.inversion_settings.bInversion_uv_loss)
            println("    inversion_bInversion_bound_loss = ", settings.inversion_settings.bInversion_bound_loss)
            println("    inversion_bInversion_slope_loss = ", settings.inversion_settings.bInversion_slope_loss)
            
            println("    inversion_bInversion_parameter_value_lower_bound = ", settings.inversion_settings.parameter_value_lower_bound)
            println("    inversion_bInversion_parameter_value_upper_bound = ", settings.inversion_settings.parameter_value_upper_bound)
            println("    inversion_bInversion_parameter_slope_upper_bound = ", settings.inversion_settings.parameter_slope_upper_bound)
              
            println("    optimizer = ", settings.inversion_settings.optimizer)
            println("    learning_rate = ", settings.inversion_settings.learning_rate)
            println("    max_iterations = ", settings.inversion_settings.max_iterations)
            println("    save_frequency = ", settings.inversion_settings.save_frequency)
            println("    save_checkpoint = ", settings.inversion_settings.save_checkpoint)
            println("    checkpoint_frequency = ", settings.inversion_settings.checkpoint_frequency)
            println("    inversion_sensealg = ", settings.inversion_settings.inversion_sensealg)
            println("    inversion_truth_file_name = ", settings.inversion_settings.inversion_truth_file_name)
            println("    inversion_forward_simulation_initial_condition_options = ", settings.inversion_settings.forward_simulation_initial_condition_options)
            println("    inversion_forward_simulation_initial_condition_file_name = ", settings.inversion_settings.forward_simulation_initial_condition_file_name)
            println("    inversion_forward_simulation_initial_condition_constant_values = ", settings.inversion_settings.forward_simulation_initial_condition_constant_values)
            println("    ode_solver = ", settings.inversion_settings.ode_solver)
            println("    ode_solver_adaptive = ", settings.inversion_settings.ode_solver_adaptive)
            println("    ode_solver_sensealg = ", settings.inversion_settings.ode_solver_sensealg)
            println("    ode_solver_b_jac_sparsity = ", settings.inversion_settings.ode_solver_b_jac_sparsity)
            println("    ode_solver_nSave = ", settings.inversion_settings.ode_solver_nSave)
            println("    save_file_name = ", settings.inversion_settings.save_file_name)
            println("    save_loss_history_file_name = ", settings.inversion_settings.save_loss_history_file_name)
            println("    save_parameters_history_file_name = ", settings.inversion_settings.save_parameters_history_file_name)
        else
            println("No inversion is to be performed.")
        end

        if settings.bPerform_Sensitivity_Analysis
            println("Sensitivity analysis is to be performed.")

            println("Sensitivity analysis options:")
            println("    active_param_names = ", settings.sensitivity_analysis_settings.active_param_names)

            if settings.sensitivity_analysis_settings.parameter_values_options == "from_file"
                println("    parameter_values_from_file = ", settings.sensitivity_analysis_settings.parameter_values_from_file)
            elseif settings.sensitivity_analysis_settings.parameter_values_options == "constant"
                println("    zb_values = ", settings.sensitivity_analysis_settings.zb_values)
                println("    ManningN_values = ", settings.sensitivity_analysis_settings.ManningN_values)
                println("    inlet_discharge_values = ", settings.sensitivity_analysis_settings.inlet_discharge_values)
            else
                error("Invalid sensitivity_parameter_values_options: $(settings.sensitivity_analysis_settings.parameter_values_options). Supported options: from_file, constant.")
            end

            println("    ode_solver = ", settings.sensitivity_analysis_settings.ode_solver)
            println("    ode_solver_adaptive = ", settings.sensitivity_analysis_settings.ode_solver_adaptive)
            println("    ode_solver_b_jac_sparsity = ", settings.sensitivity_analysis_settings.ode_solver_b_jac_sparsity)
            println("    ode_solver_nSave = ", settings.sensitivity_analysis_settings.ode_solver_nSave)
            println("    save_file_name = ", settings.sensitivity_analysis_settings.save_file_name)
        else
            println("No sensitivity analysis is to be performed.")
        end

        if settings.bPerform_UDE
            println("UDE is to be performed.")
            println("UDE options:")
            println("    UDE_mode = ", settings.UDE_settings.UDE_mode)
            println("    UDE_choice = ", settings.UDE_settings.UDE_choice)
            println("    UDE_NN_config = ", settings.UDE_settings.UDE_NN_config)
            println("    UDE_bWSE_loss = ", settings.UDE_settings.UDE_bWSE_loss)
            println("    UDE_b_uv_loss = ", settings.UDE_settings.UDE_b_uv_loss)
            println("    UDE_optimizer = ", settings.UDE_settings.UDE_optimizer)
            println("    UDE_learning_rate = ", settings.UDE_settings.UDE_learning_rate)
            println("    UDE_max_iterations = ", settings.UDE_settings.UDE_max_iterations)
            println("    UDE_training_save_frequency = ", settings.UDE_settings.UDE_training_save_frequency)
            println("    UDE_training_save_checkpoint = ", settings.UDE_settings.UDE_training_save_checkpoint)
            println("    UDE_training_checkpoint_frequency = ", settings.UDE_settings.UDE_training_checkpoint_frequency)
            println("    UDE_sensealg = ", settings.UDE_settings.UDE_sensealg)
            println("    UDE_truth_file_name = ", settings.UDE_settings.UDE_truth_file_name)
            println("    UDE_forward_simulation_initial_condition_options = ", settings.UDE_settings.UDE_forward_simulation_initial_condition_options)
            println("    UDE_forward_simulation_initial_condition_file_name = ", settings.UDE_settings.UDE_forward_simulation_initial_condition_file_name)
            println("    UDE_forward_simulation_initial_condition_constant_values = ", settings.UDE_settings.UDE_forward_simulation_initial_condition_constant_values)
            println("    ode_solver = ", settings.UDE_settings.UDE_ode_solver)
            println("    ode_solver_adaptive = ", settings.UDE_settings.UDE_ode_solver_adaptive)
            println("    ode_solver_sensealg = ", settings.UDE_settings.UDE_ode_solver_sensealg)
            println("    ode_solver_b_jac_sparsity = ", settings.UDE_settings.UDE_ode_solver_b_jac_sparsity)
            println("    ode_solver_nSave = ", settings.UDE_settings.UDE_ode_solver_nSave)
            println("    UDE_training_save_file_name = ", settings.UDE_settings.UDE_training_save_file_name)
            println("    UDE_training_save_loss_history_file_name = ", settings.UDE_settings.UDE_training_save_loss_history_file_name)
            println("    UDE_training_save_parameters_history_file_name = ", settings.UDE_settings.UDE_training_save_parameters_history_file_name)
            println("    UDE_training_save_NN_weights_state_file_name = ", settings.UDE_settings.UDE_training_save_NN_weights_state_file_name)
            println("    UDE_inference_save_file_name = ", settings.UDE_settings.UDE_inference_save_file_name)
        else
            println("No UDE is to be performed.")
        end

        println("--------------------------------")
    end

    return settings
end

# Helper functions to parse specific settings
function parse_forward_settings(options::Dict, control_file_dir::String)
    # Handle initial condition values from file if specified
    initial_condition_values = nothing
    if options["forward_simulation_initial_condition_options"] == "from_file"
        # Use the current directory to find the file
        println("Current working directory: ", pwd())
        file_path = joinpath(control_file_dir, options["forward_simulation_initial_condition_file_name"])
        println("Attempting to open file: ", file_path)

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

function parse_inversion_settings(options::Dict, control_file_dir::String)
    # Handle parameter initial values from file if specified
    parameter_values_from_file = nothing
    if options["inversion_parameter_initial_values_options"] == "from_file"
        file_path = joinpath(control_file_dir, options["inversion_parameter_initial_values_file_name"])
        parameter_values_from_file = JSON3.read(open(file_path), Dict)
    end

    forward_simulation_initial_condition_values_from_file = nothing
    if options["inversion_forward_simulation_initial_condition_options"] == "from_file"
        file_path = joinpath(control_file_dir, options["inversion_forward_simulation_initial_condition_file_name"])
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
        options["inversion_bInversion_WSE_loss"],                    # bInversion_WSE_loss
        options["inversion_bInversion_uv_loss"],                    # bInversion_uv_loss
        options["inversion_bInversion_bound_loss"],                  # bInversion_bound_loss
        options["inversion_bInversion_slope_loss"],                  # bInversion_slope_loss
        Float64.(options["inversion_parameter_value_lower_bound"]),  # parameter_value_lower_bound
        Float64.(options["inversion_parameter_value_upper_bound"]),  # parameter_value_upper_bound
        Float64.(options["inversion_parameter_slope_upper_bound"]),  # parameter_slope_upper_bound
        options["inversion_optimizer"],                        # optimizer
        Float64(options["inversion_learning_rate"]),           # learning_rate
        Int(options["inversion_max_iterations"]),              # max_iterations
        Int(options["inversion_save_frequency"]),           # save_frequency
        options["inversion_save_checkpoint"],                # save_checkpoint
        Int(options["inversion_checkpoint_frequency"]),     # checkpoint_frequency
        options["inversion_sensealg"],                       # inversion_sensealg
        options["inversion_truth_file_name"],         # inversion_truth_file_name
        options["inversion_forward_simulation_initial_condition_options"], # forward_simulation_initial_condition_options
        options["inversion_forward_simulation_initial_condition_file_name"], # forward_simulation_initial_condition_file_name
        Float64.(options["inversion_forward_simulation_initial_condition_constant_values"]), # forward_simulation_initial_condition_constant_values,
        forward_simulation_initial_condition_values_from_file, # forward_simulation_initial_condition_values_from_file
        options["inversion_ode_solver_options"]["ode_solver"], # ode_solver
        options["inversion_ode_solver_options"]["ode_solver_adaptive"], # ode_solver_adaptive
        options["inversion_ode_solver_options"]["ode_solver_sensealg"], # ode_solver_sensealg
        options["inversion_ode_solver_options"]["ode_solver_b_jac_sparsity"], # ode_solver_b_jac_sparsity
        Int(options["inversion_ode_solver_options"]["ode_solver_nSave"]), # ode_solver_nSave
        options["inversion_save_file_name"],                   # save_file_name
        get(options, "inversion_save_loss_history_file_name", "loss_history.jld2"), # save_loss_history_file_name
        get(options, "inversion_save_parameters_history_file_name", "parameters_history.csv") # save_parameters_history_file_name
    )
end

function parse_sensitivity_analysis_settings(options::Dict, control_file_dir::String)
    # Handle parameter initial values from file if specified
    parameter_values_from_file = nothing
    if options["sensitivity_parameter_values_options"] == "from_file"
        file_path = joinpath(control_file_dir, options["sensitivity_parameter_values_file_name"])
        parameter_values_from_file = JSON3.read(open(file_path), Dict)
    end

    forward_simulation_initial_condition_values_from_file = nothing
    if options["sensitivity_forward_simulation_initial_condition_options"] == "from_file"
        file_path = joinpath(control_file_dir, options["sensitivity_forward_simulation_initial_condition_file_name"])
        forward_simulation_initial_condition_values_from_file = JSON3.read(open(file_path), Dict)
    end

    SensitivityAnalysisSettings(
        String.(options["active_param_names"]),
        options["sensitivity_parameter_values_options"],
        options["sensitivity_parameter_values_file_name"],
        parameter_values_from_file,
        Float64.(options["sensitivity_zb_values"]),
        Float64.(options["sensitivity_ManningN_values"]),
        Float64.(options["sensitivity_inlet_discharge_values"]),
        options["sensitivity_forward_simulation_initial_condition_options"], # forward_simulation_initial_condition_options
        options["sensitivity_forward_simulation_initial_condition_file_name"], # forward_simulation_initial_condition_file_name
        Float64.(options["sensitivity_forward_simulation_initial_condition_constant_values"]), # forward_simulation_initial_condition_constant_values,
        forward_simulation_initial_condition_values_from_file, # forward_simulation_initial_condition_values_from_file
        options["sensitivity_ode_solver_options"]["ode_solver"],
        options["sensitivity_ode_solver_options"]["ode_solver_adaptive"],
        options["sensitivity_ode_solver_options"]["ode_solver_b_jac_sparsity"],
        Int(options["sensitivity_ode_solver_options"]["ode_solver_nSave"]),
        options["sensitivity_save_file_name"]
    )
end

function parse_UDE_settings(options::Dict, control_file_dir::String)

    #sanity check
    if !haskey(options, "UDE_mode")
        error("UDE_mode is not present in the control file.")
    end

    if options["UDE_mode"] != "training" && options["UDE_mode"] != "inference"
        error("Invalid UDE_mode: $(options["UDE_mode"]). Supported options: training, inference.")
    end

    # if it is for inference, we need to check if the NN_weights_state_file_name is provided
    if options["UDE_mode"] == "inference"
        if options["UDE_choice"]["how_to_initialize_NN"] !== "from_pretrained"
            error("For inference, how_to_initialize_NN must be from_pretrained. Does not make sense to use a random NN for inference.")
        end

        if !haskey(options["UDE_NN_config"], "NN_weights_state_file_name")
            error("NN_weights_state_file_name is not present in the control file.")
        end
    end

    if !haskey(options, "UDE_choice")
        error("UDE_choice is not present in the control file.")
    end

    #for the option of ManningN_h, the input dimension must be 1 (h) and the output dimension must be 1 (ManningN)
    if options["UDE_choice"] == "ManningN_h" 

        input_dim = options["UDE_NN_config"]["input_dim"]
        output_dim = options["UDE_NN_config"]["output_dim"]

        if input_dim != 1 || output_dim != 1
            error("For ManningN_h, the input dimension must be 1 and the output dimension must be 1.")
        end

    end

    #for the option of FlowResistance, the input dimension must be 3 (h, hu, hv) and the output dimension must be 1 (tau_b)
    if options["UDE_choice"] == "FlowResistance" 

        input_dim = options["UDE_NN_config"]["input_dim"]
        output_dim = options["UDE_NN_config"]["output_dim"]

        if input_dim != 3 || output_dim != 1
            error("For FlowResistance, the input dimension must be 3 and the output dimension must be 1.")
        end
    end

    forward_simulation_initial_condition_values_from_file = nothing
    if options["UDE_forward_simulation_initial_condition_options"] == "from_file"
        file_path = joinpath(control_file_dir, options["UDE_forward_simulation_initial_condition_file_name"])
        forward_simulation_initial_condition_values_from_file = JSON3.read(open(file_path), Dict)
    end

    UDESettings(
        String(options["UDE_mode"]),
        String(options["UDE_choice"]),
        Dict(options["UDE_NN_config"]),
        Bool(options["UDE_bWSE_loss"]),
        Bool(options["UDE_b_uv_loss"]),
        String(options["UDE_optimizer"]),
        Float64(options["UDE_learning_rate"]),
        Int(options["UDE_max_iterations"]),
        Int(options["UDE_training_save_frequency"]),
        Bool(options["UDE_training_save_checkpoint"]),
        Int(options["UDE_training_checkpoint_frequency"]),
        String(options["UDE_sensealg"]),
        String(options["UDE_truth_file_name"]),
        String(options["UDE_forward_simulation_initial_condition_options"]),
        String(options["UDE_forward_simulation_initial_condition_file_name"]),
        Float64.(options["UDE_forward_simulation_initial_condition_constant_values"]),
        forward_simulation_initial_condition_values_from_file,
        String(options["UDE_ode_solver_options"]["ode_solver"]),
        Bool(options["UDE_ode_solver_options"]["ode_solver_adaptive"]),
        String(options["UDE_ode_solver_options"]["ode_solver_sensealg"]),
        Bool(options["UDE_ode_solver_options"]["ode_solver_b_jac_sparsity"]),
        Int(options["UDE_ode_solver_options"]["ode_solver_nSave"]),
        String(options["UDE_training_save_file_name"]),
        String(options["UDE_training_save_loss_history_file_name"]),
        String(options["UDE_training_save_parameters_history_file_name"]),
        String(options["UDE_training_save_NN_weights_state_file_name"]),
        String(options["UDE_inference_save_file_name"])
    )
end


