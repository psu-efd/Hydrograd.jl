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
    adaptive::Bool
    b_jac_sparsity::Bool
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
    solution_truth_file_name::String
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
    ControlSettings(
        control_dict["bVerbose"],
        control_dict["control_variables"]["srhhydro_file_name"],
        control_dict["control_variables"]["bPerform_Forward_Simulation"],
        control_dict["control_variables"]["bPerform_Inversion"],
        control_dict["control_variables"]["bPerform_Sensitivity_Analysis"],
        time_settings,
        forward_settings,
        inversion_settings
    )
end

# Helper functions to parse specific settings
function parse_forward_settings(options::Dict)
    # Implementation of forward settings parsing
    ForwardSimulationSettings(...)
end

function parse_inversion_settings(options::Dict)
    # Implementation of inversion settings parsing
    InversionSettings(...)
end