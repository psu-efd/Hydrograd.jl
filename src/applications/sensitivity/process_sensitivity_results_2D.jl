#Functions to process the sensitivity analysis results: 
#The process_sensitivity_results function is used to plot the sensitivity analysis results.

function postprocess_sensitivity_results_swe_2D(settings, my_mesh_2D, nodeCoordinates, params_vector, active_param_name, case_path)

    #get the number of cells
    n_cells = my_mesh_2D.numOfCells

    #get the length of the active parameter
    nParams = length(params_vector)

    #load the sensitivity analysis results
    sensitivity_results_file_name = joinpath(case_path, settings.sensitivity_analysis_settings.save_file_name)
    sensitivity_results = load(sensitivity_results_file_name)

    #extract the sensitivity analysis results
    sensitivity = sensitivity_results["sensitivity"]

    #if settings.bVerbose
        @show size(sensitivity)
    #end
    
    sensitivity_size = size(sensitivity)  #(n_cells * 3, n_params); each cell has 3 solution variables: h, hu, hv. Each cell has 3 rows in the sol array.

    @assert sensitivity_size[1] == n_cells * 3
    @assert sensitivity_size[2] == nParams    

    # loop over each parameter and save the sensitivity results to vtk
    for i in 1:nParams

        #initialize the sensitivity arrays
        dh_dparam = zeros(n_cells)
        dhu_dparam = zeros(n_cells)
        dhv_dparam = zeros(n_cells)

        #loop over each cell
        for j in 1:n_cells
            dh_dparam[j]  = sensitivity[j, i]     # h sensitivity
            dhu_dparam[j] = sensitivity[n_cells+j, i]   # hu sensitivity
            dhv_dparam[j] = sensitivity[2*n_cells+j, i]   # hv sensitivity
        end

        field_name = "parameter_number"
        field_type = "integer"
        field_value = i
   
        vector_data = []
        vector_names = []

        scalar_data = [dh_dparam, dhu_dparam, dhv_dparam]
        scalar_names = ["dh_dparam", "dhu_dparam", "dhv_dparam"]

        file_path = joinpath(case_path, "sensitivity_results_$(active_param_name)_$i.vtk")

        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, 
            field_name, field_type, field_value, 
            scalar_data, scalar_names, vector_data, vector_names)    

        println("       Sensitivity results for parameter $(active_param_name) = $i) is saved to ", file_path)

        #save the sensitivity results in a JSON file
        open(joinpath(case_path, "sensitivity_results_$(active_param_name)_$i.json"), "w") do io
            JSON3.pretty(io, Dict("dh_dparam" => dh_dparam, 
                                  "dhu_dparam" => dhu_dparam,
                                  "dhv_dparam" => dhv_dparam,
                                  "parameter_name" => string(active_param_name, "_", i),
                                  "parameter_value" => params_vector[i]))
            println(io)
        end
    end

end


