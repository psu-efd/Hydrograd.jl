
#preprocess: create a ComponentArray for the model parameters for 2D shallow water equations
#The complete model parameters include:
#1. Manning's n
#2. bed elevation
#3. Inlet discharge (constant for now)
#4. xxx (for future use)
function preprocess_model_parameters_2D(settings, zb_cells_param, ManningN_list_param, inlet_discharges_param)

    # Create ComponentArray for parameters
    params_array = ComponentArray(zb_cells_param=zb_cells_param, ManningN_list_param=ManningN_list_param, inlet_discharges_param=inlet_discharges_param)
    active_params = []

    #for forward simulation, all parameters are active (not in the sense of inversion, but for the sake of passing parameters to the forward simulation)
    #the parameter values are from SRH-2D data
    if settings.bPerform_Forward_Simulation   
        active_params = [:zb_cells_param, :ManningN_list_param, :inlet_discharges_param]
    end

    if settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis
        if length(settings.inversion_settings.active_param_names) != 1
            error("Currently only support one active parameter for now. Active parameter names: $active_param_names. Supported active parameter names: zb, ManningN, Q")
        end

        #get the active parameter name (only one active parameter is supported for now)
        active_param_name = settings.inversion_settings.active_param_names[1]

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
