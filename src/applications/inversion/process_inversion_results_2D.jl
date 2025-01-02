#Functions to process the inversion results: 
#The process_inversion_results function is used to process the inversion results.

function  postprocess_inversion_results_swe_2D(settings, my_mesh_2D, nodeCoordinates, zb_cell_truth, h_truth, u_truth, v_truth, WSE_truth, case_path)

    inversion_results_file_name = joinpath(case_path, settings.inversion_settings.save_file_name)

    inversion_save_loss_history_file_name = settings.inversion_settings.save_loss_history_file_name


    inversion_results = load(inversion_results_file_name)
    
    LOSS = inversion_results["LOSS"]
    PRED = inversion_results["PRED"]
    PARS = inversion_results["PARS"]
    
    #extract the Losses
    iter_numbers = collect(range(1,size(LOSS)[1]))

    loss_total = []
    loss_pred = []
    loss_pred_WSE = []
    loss_pred_u = []
    loss_slope = []

    for curLoss in LOSS
        append!(loss_total, curLoss[1])
        append!(loss_pred, curLoss[2])
        append!(loss_pred_WSE, curLoss[3])
        append!(loss_pred_u, curLoss[4])
        append!(loss_slope, curLoss[5])
    end 

    loss_history_df = DataFrame(iter_numbers = iter_numbers, 
                                loss_total = loss_total, 
                                loss_pred = loss_pred, 
                                loss_pred_WSE = loss_pred_WSE, 
                                loss_pred_u = loss_pred_u, 
                                loss_slope = loss_slope)

    #save loss history to a file
    CSV.write(joinpath(case_path, inversion_save_loss_history_file_name), loss_history_df)

    #save the results of each iteration to vtk 
    for (i, (curPred, curPars)) in enumerate(zip(PRED, PARS))
        h_i = curPred
        
        #if the inversion is for zb, then curPars is the zb_cells_param. Otherwise, curPar is something else such as the Manning's n or the inlet discharges. 
        if settings.inversion_settings.active_param_names == ["zb"]
            zb_i = curPars
        else
            zb_i = zb_cell_truth
        end

        WSE_i = h_i .+ zb_i

        vector_data = [] 
        vector_names = []

        scalar_data = [h_i, WSE_i, zb_i, h_truth, WSE_truth, zb_cell_truth]
        scalar_names = ["h", "WSE", "zb", "h_truth", "WSE_truth", "zb_truth"]
        
        file_path = joinpath(case_path, "inversion_results_iter_$i.vtk" ) 
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, 
                         scalar_data, scalar_names, vector_data, vector_names)    
    
        println("inversion_results_iter_$i.vtk is saved to ", file_path)
    end

    
end


