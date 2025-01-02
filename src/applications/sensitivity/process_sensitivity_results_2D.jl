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
    sol = sensitivity_results["sol"]

    if settings.bVerbose
        @show size(sol)
    end
    
    sol_size = size(sol)  #(n_cells * 3, n_params); each cell has 3 solution variables: h, hu, hv. Each cell has 3 rows in the sol array.

    @assert sol_size[1] == n_cells * 3
    @assert sol_size[2] == nParams    

    # loop over each parameter and save the sensitivity results to vtk
    for i in 1:nParams

        #initialize the sensitivity arrays
        dh_dparam = zeros(n_cells)
        dhu_dparam = zeros(n_cells)
        dhv_dparam = zeros(n_cells)

        #loop over each cell
        for j in 1:n_cells
            idx = (j - 1) * 3 + 1           # Starting row index for each cell
            dh_dparam[j]  = sol[idx, i]     # h sensitivity
            dhu_dparam[j] = sol[idx+1, i]   # hu sensitivity
            dhv_dparam[j] = sol[idx+2, i]   # hv sensitivity
        end
   
        vector_data = []
        vector_names = []

        scalar_data = [dh_dparam, dhu_dparam, dhv_dparam]
        scalar_names = ["dh_dparam", "dhu_dparam", "dhv_dparam"]

        file_path = joinpath(case_path, "sensitivity_results_iParam_$i.vtk")
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount,
            scalar_data, scalar_names, vector_data, vector_names)

        println("       Sensitivity results for parameter $i is saved to ", file_path)
    end

end


