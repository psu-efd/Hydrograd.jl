#Functions to process the forward simulation results: 
#The process_forward_simulation_results function is used to plot the forward simulation results in Python.

function  postprocess_forward_simulation_results_swe_2D(settings, my_mesh_2D, zb_cell_truth, 
    ManningN_cell_truth, inlet_discharges_truth, zb_cells, nodeCoordinates, case_path)

    forward_simulation_results_file_name = joinpath(case_path, settings.forward_settings.save_file_name)
    forward_simulation_results = load(forward_simulation_results_file_name)["sol_data"]["u"]  #get the solution data (state variable)

    #@show typeof(forward_simulation_results)
    #@show size(forward_simulation_results)
    #@show size(forward_simulation_results[end][:,:])
    #@show forward_simulation_results[end][:,:]

    #save the simulation results (h, u, v) at the last time step to a json file (to be used as ground truth for inversion)
    h_truth = vec(forward_simulation_results)[end][:, 1]
    wse_truth = h_truth .+ zb_cell_truth
    u_truth = vec(forward_simulation_results)[end][:, 2] ./ vec(forward_simulation_results)[end][:, 1]
    v_truth = vec(forward_simulation_results)[end][:, 3] ./ vec(forward_simulation_results)[end][:, 1]
    open(joinpath(case_path, settings.forward_settings.save_solution_truth_file_name), "w") do io
        JSON3.pretty(io, Dict("wse_truth" => wse_truth, "h_truth" => h_truth, "u_truth" => u_truth, "v_truth" => v_truth, "zb_cell_truth" => zb_cell_truth, 
                            "ManningN_cell_truth" => ManningN_cell_truth, "inlet_discharges_truth" => inlet_discharges_truth))
        println(io)
    end

    swe_2D_save_results_SciML(forward_simulation_results, my_mesh_2D, nodeCoordinates, zb_cells, case_path)
    
end


