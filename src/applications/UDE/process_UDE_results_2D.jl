#Functions to process the UDE results: 
#The process_UDE_results function is used to process the UDE results.

function  postprocess_UDE_training_results_swe_2D(swe2d_extra_params, zb_cell_truth, h_truth, u_truth, v_truth, WSE_truth, 
                                 ManningN_cells_truth, friction_x_truth, friction_y_truth)

    settings = swe2d_extra_params.settings
    my_mesh_2D = swe2d_extra_params.my_mesh_2D
    nodeCoordinates = swe2d_extra_params.nodeCoordinates
    case_path = swe2d_extra_params.case_path
    ManningN_cells = swe2d_extra_params.ManningN_cells

    ude_model = swe2d_extra_params.ude_model
    ude_model_state = swe2d_extra_params.ude_model_state

    #friction magnitudes
    friction_magnitudes_truth = sqrt.(friction_x_truth.^2 .+ friction_y_truth.^2)

    UDE_training_results_file_name = joinpath(case_path, settings.UDE_settings.UDE_training_save_file_name)
    UDE_training_save_loss_history_file_name = settings.UDE_settings.UDE_training_save_loss_history_file_name

    UDE_training_results = load(UDE_training_results_file_name)
    
    ITER = UDE_training_results["ITER"]
    LOSS = UDE_training_results["LOSS"]
    PARS = UDE_training_results["PARS"]

    #println("ITER: ", ITER)
    #println("LOSS: ", LOSS)
    #println("PARS: ", PARS)

    #Load the UDE model, parameters (only the last iteration), and state
    ude_model_params = UDE_training_results["ude_model_params"]
    ude_model_state = UDE_training_results["ude_model_state"]

    #@show typeof(ude_model_params)
    #@show ude_model_params
    
    #extract the Losses
    loss_total = []
    loss_pred_WSE = []
    loss_pred_uv = []

    for curLoss in LOSS
        append!(loss_total, curLoss[1])
        append!(loss_pred_WSE, curLoss[2])        
        append!(loss_pred_uv, curLoss[3])
    end 

    loss_history_df = DataFrame(iter_numbers = ITER, 
                                loss_total = loss_total, 
                                loss_pred_WSE = loss_pred_WSE, 
                                loss_pred_uv = loss_pred_uv)

    #save loss history to a file
    CSV.write(joinpath(case_path, UDE_training_save_loss_history_file_name), loss_history_df)

    #@show typeof(PARS)
    #@show size(PARS)
    #@show PARS

    #save the results of each iteration to vtk 
    for (i, (curPred, curPars)) in enumerate(zip(PRED, PARS))
        println("Processing UDE training results for iteration $(ITER[i])")

        #@show typeof(curPars)
        #@show curPars
        

        #If training the UDE and the UDE choice is ManningN_h
        if settings.UDE_settings.UDE_choice == "ManningN_h"
            h_i = curPred[:, 1]
            # Reconstruct the parameters dictionary for the neural network
            #reconstructed_params = ComponentArrays.structuralize(curPars, swe2d_extra_params.ude_model_params)
            ManningN_cells = update_ManningN_UDE(h_i, ude_model, curPars, ude_model_state, Float64.(settings.UDE_settings.UDE_NN_config["h_bounds"]), my_mesh_2D.numOfCells)
        end

        #If training the UDE and the UDE choice is FlowResistance
        if settings.UDE_settings.UDE_choice == "FlowResistance"
            h_i = curPred[:, 1]
            q_x_i = curPred[:, 2]
            q_y_i = curPred[:, 3]
            friction_magnitudes = update_FlowResistance_UDE(h_i, q_x_i, q_y_i, ude_model, curPars, ude_model_state, Float64.(settings.UDE_settings.UDE_NN_config["h_bounds"]), 
                       Float64.(settings.UDE_settings.UDE_NN_config["velocity_magnitude_bounds"]), my_mesh_2D.numOfCells)
        end

        WSE_i = h_i .+ zb_cell_truth

        field_name = "iter_number"
        field_type = "integer"
        field_value = i

        vector_data = [] 
        vector_names = []

        if settings.UDE_settings.UDE_choice == "ManningN_h"
            scalar_data = [h_i, WSE_i, h_truth, WSE_truth, zb_cell_truth, ManningN_cells, ManningN_cells_truth]
            scalar_names = ["h", "WSE", "h_truth", "WSE_truth", "zb_truth", "ManningN", "ManningN_truth"]
        elseif settings.UDE_settings.UDE_choice == "FlowResistance"
            scalar_data = [h_i, WSE_i, h_truth, WSE_truth, zb_cell_truth, friction_magnitudes, friction_magnitudes_truth]
            scalar_names = ["h", "WSE", "h_truth", "WSE_truth", "zb_truth", "friction_magnitudes", "friction_magnitudes_truth"]
        else
            error("UDE choice $(settings.UDE_settings.UDE_choice) is not valid. Please check the UDE choice in the control file.")
        end

        file_path = joinpath(case_path, "UDE_training_results_iter_$i.vtk" ) 
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, 
                         field_name, field_type, field_value, scalar_data, scalar_names, vector_data, vector_names)    
    
        println("       UDE training results for iteration $i is saved to ", file_path)
    end

    
end


