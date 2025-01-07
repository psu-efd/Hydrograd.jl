#Functions to process the forward simulation results: 
#The process_forward_simulation_results function is used to process the forward simulation results.

function  postprocess_forward_simulation_results_swe_2D(swe2d_extra_params, zb_cell_truth, ManningN_zone_values_truth, inlet_discharges_truth)

    settings = swe2d_extra_params.settings
    my_mesh_2D = swe2d_extra_params.my_mesh_2D
    nodeCoordinates = swe2d_extra_params.nodeCoordinates
    zb_cells = swe2d_extra_params.zb_cells
    case_path = swe2d_extra_params.case_path

    ManningN_cells = swe2d_extra_params.ManningN_cells

    forward_simulation_results_file_name = joinpath(case_path, settings.forward_settings.save_file_name)
    forward_simulation_results = load(forward_simulation_results_file_name)["sol_data"]["u"]  #get the solution data (state variable)

    #save the simulation results (h, u, v) at the last time step to a json file (to be used as ground truth for inversion)
    h_truth = vec(forward_simulation_results)[end][:, 1]
    wse_truth = h_truth .+ zb_cell_truth
    q_x_truth = vec(forward_simulation_results)[end][:, 2]
    q_y_truth = vec(forward_simulation_results)[end][:, 3]
    u_truth = q_x_truth ./ h_truth
    v_truth = q_y_truth ./ h_truth

    #If ManningN_option is variable_as_function_of_h, update ManningN_cell based on the ManningN_function_type and ManningN_function_parameters
    #if settings.forward_settings.ManningN_option == "variable_as_function_of_h"
    #    ManningN_cells = update_ManningN_forward_simulation(h_truth, settings)
    #end

    #compute the friction terms
    friction_x_truth, friction_y_truth = compute_friction_terms(settings, h_truth, q_x_truth, q_y_truth, swe2d_extra_params.ManningN_cells, 
                 nothing, swe2d_extra_params, my_mesh_2D, swe2d_extra_params.swe_2D_constants.g, swe2d_extra_params.swe_2D_constants.h_small)

    #save the simulation results as truth to a json file
    open(joinpath(case_path, settings.forward_settings.save_solution_truth_file_name), "w") do io
        JSON3.pretty(io, Dict("wse_truth" => wse_truth, "h_truth" => h_truth, "u_truth" => u_truth, "v_truth" => v_truth, "zb_cell_truth" => zb_cell_truth, 
                            "ManningN_cells_truth" => ManningN_cells, 
                            "friction_x_truth" => friction_x_truth, "friction_y_truth" => friction_y_truth,
                            "ManningN_zone_values_truth" => ManningN_zone_values_truth, "inlet_discharges_truth" => inlet_discharges_truth))
        println(io)
    end

    #save the simulation results to vtk files
    swe_2D_save_results_SciML(swe2d_extra_params, forward_simulation_results, friction_x_truth, friction_y_truth)
    
end


