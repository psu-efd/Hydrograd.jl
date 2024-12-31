#Functions to process the forward simulation results: 
#The process_forward_simulation_results function is used to plot the forward simulation results in Python.

function  postprocess_forward_simulation_results_swe_2D(settings, my_mesh_2D, total_water_volume, zb_cell_truth, 
    ManningN_cell_truth, inlet_discharges_truth, zb_cells, nodeCoordinates, case_path)

    forward_simulation_results_file_name = joinpath(case_path, settings.forward_settings.save_file_name)
    forward_simulation_results = load(forward_simulation_results_file_name)["sol"]

    #save the simulation results (h, u, v) at the last time step to a json file (to be used as ground truth for inversion)
    h_truth = Array(forward_simulation_results)[:, 1, end]
    wse_truth = h_truth .+ zb_cell_truth
    u_truth = Array(forward_simulation_results)[:, 2, end] ./ Array(forward_simulation_results)[:, 1, end]
    v_truth = Array(forward_simulation_results)[:, 3, end] ./ Array(forward_simulation_results)[:, 1, end]
    open(joinpath(case_path, settings.forward_settings.save_solution_truth_file_name), "w") do io
        JSON3.pretty(io, Dict("wse_truth" => wse_truth, "h_truth" => h_truth, "u_truth" => u_truth, "v_truth" => v_truth, "zb_cell_truth" => zb_cell_truth, 
                            "ManningN_cell_truth" => ManningN_cell_truth, "inlet_discharges_truth" => inlet_discharges_truth))
        println(io)
    end

    swe_2D_save_results_SciML(forward_simulation_results, total_water_volume, my_mesh_2D, nodeCoordinates, zb_cells, case_path)
    
end


