
#preprocess: create a ComponentArray for the model parameters for 2D shallow water equations
#The complete model parameters include:
#1. Manning's n
#2. bed elevation
#3. Inlet discharge (constant for now)
#4. xxx (for future use)
function preprocess_model_parameters_2D(bPerform_Forward_Simulation, bPerform_Inversion, bPerform_Sensitivity_Analysis, zb_cells_param, 
    ManningN_list_param, inlet_discharges_param, active_param_names)

    # Create ComponentArray for parameters
    params_array = ComponentArray(zb_cells_param=zb_cells_param, ManningN_list_param=ManningN_list_param, inlet_discharges_param=inlet_discharges_param)
    active_params = []

    if bPerform_Inversion || bPerform_Sensitivity_Analysis
        if length(active_param_names) != 1
            error("Currently only support one active parameter for now. Active parameter names: $active_param_names. Supported active parameter names: zb, ManningN, Q")
        end

        #get the active parameter name
        active_param_name = active_param_names[1]

        if active_param_name == "zb"
            active_params = [:zb_cells_param]
        elseif active_param_name == "ManningN"
            active_params = [:ManningN_list_param]
        elseif active_param_name == "Q"
            active_params = [:inlet_discharges_param]
        else
            error("Invalid active parameter name: $active_param_name")
        end
    end

    return params_array, active_params
end
