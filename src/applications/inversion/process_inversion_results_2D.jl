#Functions to process the inversion results: 
#The process_inversion_results function is used to process the inversion results.

function  postprocess_inversion_results_swe_2D(settings, my_mesh_2D, nodeCoordinates, zb_cell_truth, h_truth, u_truth, v_truth, WSE_truth, case_path)

    inversion_results_file_name = joinpath(case_path, settings.inversion_settings.save_file_name)
    inversion_save_loss_history_file_name = settings.inversion_settings.save_loss_history_file_name
    inversion_save_parameters_history_file_name = settings.inversion_settings.save_parameters_history_file_name


    inversion_results = load(inversion_results_file_name)
    
    ITER = inversion_results["ITER"]
    LOSS = inversion_results["LOSS"]
    #PRED = inversion_results["PRED"]
    PARS = inversion_results["PARS"]
    
    #extract the Losses
    loss_total = []
    loss_pred = []
    loss_pred_WSE = []
    loss_pred_uv = []
    loss_bound = []
    loss_slope = []

    for curLoss in LOSS
        append!(loss_total, curLoss[1])
        append!(loss_pred, curLoss[2])
        append!(loss_pred_WSE, curLoss[3])
        append!(loss_pred_uv, curLoss[4])
        append!(loss_bound, curLoss[5])
        append!(loss_slope, curLoss[6])
    end 

    loss_history_df = DataFrame(iter_numbers = ITER, 
                                loss_total = loss_total, 
                                loss_pred = loss_pred, 
                                loss_pred_WSE = loss_pred_WSE, 
                                loss_pred_uv = loss_pred_uv, 
                                loss_bound = loss_bound, 
                                loss_slope = loss_slope)

    #save loss history to a file
    CSV.write(joinpath(case_path, inversion_save_loss_history_file_name), loss_history_df)

    #@show typeof(PARS)
    #@show size(PARS)
    #@show PARS

    #save parameters to a file
    # Create column names based on parameter count
    param_name = settings.inversion_settings.active_param_names[1]
    num_params = length(PARS[1])  # Get number of parameters from first vector
    param_names = ["$(param_name)_$i" for i in 1:num_params]

    pars_df = hcat(DataFrame(iter_numbers = ITER), DataFrame(hcat(PARS...)', param_names))

    # Write to CSV
    CSV.write(joinpath(case_path, inversion_save_parameters_history_file_name), pars_df)

    #save the results of each iteration to vtk 
    for (i, (curPred, curPars)) in enumerate(zip(PRED, PARS))
        @show size(curPred)

        h_i = curPred[1:my_mesh_2D.numOfCells]
        
        #if the inversion is for zb, then curPars is the zb_cells_param. Otherwise, curPar is something else such as the Manning's n or the inlet discharges. 
        if settings.inversion_settings.active_param_names == ["zb"]
            zb_i = curPars

            #compute S0
            _, _, S0_i = interploate_zb_from_cell_to_face_and_compute_S0(my_mesh_2D, curPars)
        else
            zb_i = zb_cell_truth
        end

        WSE_i = h_i .+ zb_i

        field_name = "iter_number"
        field_type = "integer"
        field_value = i

        vector_data = [] 
        vector_names = []

        if settings.inversion_settings.active_param_names == ["zb"]
            vector_data = [S0_i]
            vector_names = ["S0"]

            scalar_data = [h_i, WSE_i, zb_i, h_truth, WSE_truth, zb_cell_truth]
            scalar_names = ["h", "WSE", "zb", "h_truth", "WSE_truth", "zb_truth"]
        else
            scalar_data = [h_i, WSE_i, h_truth, WSE_truth]
            scalar_names = ["h", "WSE", "h_truth", "WSE_truth"]
        end
        
        file_path = joinpath(case_path, "inversion_results_iter_$i.vtk" ) 
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, 
                         field_name, field_type, field_value, scalar_data, scalar_names, vector_data, vector_names)    
    
        println("       Inversion results for iteration $i is saved to ", file_path)
    end

    
end


