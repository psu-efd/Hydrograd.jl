
#preprocess: create a ComponentArray for the model parameters for 2D shallow water equations
#The complete model parameters include:
#1. Manning's n
#2. bed elevation
#3. Inlet discharge (constant for now)
#4. xxx (for future use)
function preprocess_model_parameters_2D(zb_cells_param, ManningN_list_param, inlet_discharges_param, active_param_names)

    #currently only support one active parameter for now
    if length(active_param_names) != 1
        error("Currently only support one active parameter for now. Active parameter names: $active_param_names. Supported active parameter names: zb, ManningN, Q")
    end

    #get the active parameter name
    active_param_name = active_param_names[1]

    # Create ComponentArray for parameters
    params_array = ComponentArray(zb_cells_param=zb_cells_param, ManningN_list_param=ManningN_list_param, inlet_discharges_param=inlet_discharges_param)

    if active_param_name == "zb"    
        active_params = [:zb_cells_param]
    elseif active_param_name == "ManningN"
        active_params = [:ManningN_list_param]
    elseif active_param_name == "Q"
        active_params = [:inlet_discharges_param]
    else
        error("Invalid active parameter name: $active_param_name")
    end

    return params_array, active_params

end
