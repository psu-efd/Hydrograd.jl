
#preprocess: create a parameter vector for 2D shallow water equations
#The complete model parameter vector include:
#1. bed elevation (ncells)
#2. Manning's n (nMaterials)
#3. Inlet discharge (nInlet)
#4. xxx (for future use)
function preprocess_model_parameters_2D(settings, zb_cells_param, ManningN_list_param, inlet_discharges_param)

    if settings.bVerbose
        println("preprocess_model_parameters_2D")
    end

    # Concatenate all parameters into a single vector
    if !isnothing(inlet_discharges_param)
        params_array = vcat(zb_cells_param, ManningN_list_param, inlet_discharges_param)
    else
        params_array = vcat(zb_cells_param, ManningN_list_param)
    end

    # First calculate all the lengths and indices
    n_zb = length(zb_cells_param)
    n_manning = length(ManningN_list_param)

    if !isnothing(inlet_discharges_param)
        n_inlet = length(inlet_discharges_param)
    else
        n_inlet = 0
    end
    
    # Calculate start and end indices
    zb_start = 1
    zb_end = n_zb
    manning_start = zb_end + 1
    manning_end = manning_start + n_manning - 1

    if !isnothing(inlet_discharges_param)
        inletQ_start = manning_end + 1
        inletQ_end = inletQ_start + n_inlet - 1
    else
        inletQ_start = nothing
        inletQ_end = nothing
    end

    # Store the range for each parameter for later parameter extraction
    param_ranges = (
        zb_start = zb_start,
        zb_end = zb_end,
        manning_start = manning_start,
        manning_end = manning_end,
        inletQ_start = inletQ_start,
        inletQ_end = inletQ_end
    )

    # Determine active ranges of parameters
    active_range= if settings.bPerform_Forward_Simulation
        #for forward simulation, all parameters are active (not in the sense of inversion, but for the sake of passing parameters to the forward simulation)
        #the parameter values are from SRH-2D data
        1:length(params_array)  # All parameters active
    
    elseif settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis
        
        if length(settings.inversion_settings.active_param_names) != 1
            error("Currently only support one active parameter for now. Active parameter names: $active_param_names. Supported active parameter names: zb, ManningN, Q")
        end

        #get the active parameter name (only one active parameter is supported for now)
        active_param = settings.inversion_settings.active_param_names[1]
        if active_param == "zb"
            param_ranges.zb_start:param_ranges.zb_end
        elseif active_param == "ManningN"
            param_ranges.manning_start:param_ranges.manning_end
        elseif active_param == "Q"
            if !isnothing(inlet_discharges_param)
                param_ranges.inletQ_start:param_ranges.inletQ_end
            else
                error("No inletQ_BCs are defined in the SRH-2D data. Make no sense to perform inversion or sensitivity analysis on Q. Please check case settings.")
            end
        else
            error("Invalid active parameter name: $active_param")
        end
    end

    if settings.bVerbose
        Zygote.ignore() do
            println("params_array = ", params_array)
            println("active_range = ", active_range)
            println("param_ranges = ", param_ranges)
        end
    end

    return params_array, active_range, param_ranges
end
