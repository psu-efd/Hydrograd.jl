#Setup the parameter arrays: zb_cells_param, ManningN_list_param, inlet_discharges_param
function setup_model_parameters_2D(settings, my_mesh_2D, srh_all_Dict, zb_cells_truth, ManningN_values_truth, inlet_discharges_truth)

    #initialize the parameter arrays
    zb_cells_param = zeros(eltype(zb_cells_truth), length(zb_cells_truth))
    ManningN_list_param = zeros(eltype(ManningN_values_truth), length(ManningN_values_truth))
    inlet_discharges_param = zeros(eltype(inlet_discharges_truth), length(inlet_discharges_truth))

    active_param_name = "none"

    #For forward simulation, the parameter values are from SRH-2D data. These parameter values are not used in the forward simulation.
    #For inversion and sensitivity analysis, the parameter values are from the initial values or from a file.
    if settings.bPerform_Forward_Simulation
        zb_cells_param = deepcopy(zb_cells_truth)
        ManningN_list_param = deepcopy(ManningN_values_truth)

        if !isnothing(inlet_discharges_truth)
            inlet_discharges_param = deepcopy(inlet_discharges_truth)
        end
    elseif settings.bPerform_Inversion 
        if settings.inversion_settings.parameter_initial_values_options == "constant"

            if settings.inversion_settings.active_param_names == ["zb"]
                active_param_name = "zb"
                #bed elevation
                zb_cells_param = fill(settings.inversion_settings.zb_initial_values[1], my_mesh_2D.numOfCells)
            elseif settings.inversion_settings.active_param_names == ["ManningN"]
                active_param_name = "ManningN"
                #Manning's n
                ManningN_list_param = deepcopy(settings.inversion_settings.ManningN_initial_values)
            elseif settings.inversion_settings.active_param_names == ["Q"]
                active_param_name = "Q"
                #inlet discharges
                inlet_discharges_param = deepcopy(settings.inversion_settings.inlet_discharge_initial_values)
            else
                error("Invalid active parameter name: $(settings.inversion_settings.active_param_names). Supported active parameter names: zb, ManningN, Q")
            end
        elseif settings.inversion_settings.parameter_initial_values_options == "from_file"
            if settings.inversion_settings.active_param_names == ["zb"]
                active_param_name = "zb"
                zb_cells_param = settings.inversion_settings.parameter_initial_values_from_file["zb_cells_param"]
            elseif settings.inversion_settings.active_param_names == ["ManningN"]
                active_param_name = "ManningN"
                ManningN_list_param = settings.inversion_settings.parameter_initial_values_from_file["ManningN_list_param"]
            elseif settings.inversion_settings.active_param_names == ["Q"]
                active_param_name = "Q"
                if !isnothing(inlet_discharges_truth)
                    inlet_discharges_param = settings.inversion_settings.parameter_initial_values_from_file["inlet_discharges_param"]
                end
            else
                error("Invalid active parameter name: $(settings.inversion_settings.active_param_names). Supported active parameter names: zb, ManningN, Q")
            end
        else
            error("Invalid zb_initial_values_options: $(settings.inversion_settings.parameter_initial_values_options). Supported options: zero, from_file.")
        end

        #consistency check
        #make sure the length of zb_cells_param is the same as the number of cells
        if length(zb_cells_param) != my_mesh_2D.numOfCells
            error("The length of zb_cells_param is not the same as the number of cells.")
        end

        #make sure the length of ManningN_list_param is the same as the number of materials
        srhmat_numOfMaterials = srh_all_Dict["srhmat_numOfMaterials"]
        if length(ManningN_list_param) != srhmat_numOfMaterials
        #    error("The length of ManningN_list_param is not the same as the number of materials.")
        end

        #make sure the length of inlet_discharges_param is the same as the number of inletQ_BCs
        nInletQ_BCs = srh_all_Dict["nInletQ_BCs"]

        if !isnothing(inlet_discharges_truth)
            if length(inlet_discharges_param) != nInletQ_BCs
                error("Length mismatch: inlet_discharges_param ($(length(inlet_discharges_param))) != nInletQ_BCs ($nInletQ_BCs)")
            end
        else    #if inlet_discharges_truth is nothing, then nInletQ_BCs should be 0
            if nInletQ_BCs > 0
                error("No inletQ_BCs are defined in the SRH-2D data, but nInletQ_BCs is greater than 0. Please check the SRH-2D data.")
            end
        end
    elseif settings.bPerform_Sensitivity_Analysis
        if settings.sensitivity_analysis_settings.parameter_values_options == "constant"
            if settings.sensitivity_analysis_settings.active_param_names == ["zb"]
                active_param_name = "zb"
                #bed elevation
                zb_cells_param = fill(settings.sensitivity_analysis_settings.zb_values[1], my_mesh_2D.numOfCells)
            elseif settings.sensitivity_analysis_settings.active_param_names == ["ManningN"]
                active_param_name = "ManningN"
                #Manning's n
                ManningN_list_param = deepcopy(settings.sensitivity_analysis_settings.ManningN_values)
            elseif settings.sensitivity_analysis_settings.active_param_names == ["Q"]
                active_param_name = "Q"
                #inlet discharges
                if !isnothing(inlet_discharges_truth)
                    inlet_discharges_param = deepcopy(settings.sensitivity_analysis_settings.inlet_discharge_values)
                end
            else
                error("Invalid active parameter name: $(settings.sensitivity_analysis_settings.active_param_names). Supported active parameter names: zb, ManningN, Q")
            end
        elseif settings.sensitivity_analysis_settings.parameter_values_options == "from_file"
            if settings.sensitivity_analysis_settings.active_param_names == ["zb"]
                active_param_name = "zb"
                zb_cells_param = settings.sensitivity_analysis_settings.parameter_initial_values_from_file["zb_cells_param"]
            elseif settings.sensitivity_analysis_settings.active_param_names == ["ManningN"]
                active_param_name = "ManningN"
                ManningN_list_param = settings.sensitivity_analysis_settings.parameter_initial_values_from_file["ManningN_list_param"]
            elseif settings.sensitivity_analysis_settings.active_param_names == ["Q"]
                active_param_name = "Q"
                if !isnothing(inlet_discharges_truth)
                    inlet_discharges_param = settings.sensitivity_analysis_settings.parameter_initial_values_from_file["inlet_discharges_param"]
                end
            else
                error("Invalid active parameter name: $(settings.sensitivity_analysis_settings.active_param_names). Supported active parameter names: zb, ManningN, Q")
            end
        else
            error("Invalid parameter values option: $(settings.sensitivity_analysis_settings.parameter_values_options). Supported options: constant, from_file.")
        end

        #consistency check
        #make sure the length of zb_cells_param is the same as the number of cells
        if length(zb_cells_param) != my_mesh_2D.numOfCells
            error("The length of zb_cells_param is not the same as the number of cells.")
        end

        #make sure the length of ManningN_list_param is the same as the number of materials
        srhmat_numOfMaterials = srh_all_Dict["srhmat_numOfMaterials"]
        if length(ManningN_list_param) != srhmat_numOfMaterials
            error("The length of ManningN_list_param is not the same as the number of materials.")
        end

        #make sure the length of inlet_discharges_param is the same as the number of inletQ_BCs
        nInletQ_BCs = srh_all_Dict["nInletQ_BCs"]

        if !isnothing(inlet_discharges_truth)
            if length(inlet_discharges_param) != nInletQ_BCs
                error("Length mismatch: inlet_discharges_param ($(length(inlet_discharges_param))) != nInletQ_BCs ($nInletQ_BCs)")
            end
        else    #if inlet_discharges_truth is nothing, then nInletQ_BCs should be 0
            if nInletQ_BCs > 0
                error("No inletQ_BCs are defined in the SRH-2D data, but nInletQ_BCs is greater than 0. Please check the SRH-2D data.")
            end
        end
    end

    if active_param_name == "zb"
        params_vector = zb_cells_param
    elseif active_param_name == "ManningN"
        params_vector = ManningN_list_param
    elseif active_param_name == "Q"
        params_vector = inlet_discharges_param
    else
        params_vector = [0.0]   #not used
    end

    return params_vector, active_param_name

end