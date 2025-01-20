#Functions to process the forward simulation results: 
#The process_forward_simulation_results function is used to process the forward simulation results.

function  postprocess_forward_simulation_results_swe_2D(swe2d_extra_params::SWE2D_Extra_Parameters, 
     zb_cell_truth::Vector{T}, ManningN_zone_values_truth::Vector{T}, inlet_discharges_truth::Vector{T}) where {T<:Real}

    settings = swe2d_extra_params.settings
    my_mesh_2D = swe2d_extra_params.my_mesh_2D
    nodeCoordinates = swe2d_extra_params.nodeCoordinates
    zb_cells = swe2d_extra_params.zb_cells
    case_path = swe2d_extra_params.case_path

    ManningN_cells = swe2d_extra_params.ManningN_cells

    wstill = swe2d_extra_params.wstill    
    hstill = swe2d_extra_params.hstill

    forward_simulation_results_file_name = joinpath(case_path, settings.forward_settings.save_file_name)
    forward_simulation_results = load(forward_simulation_results_file_name)["sol_data"]["u"]  #get the solution data (state variable)

    #save the simulation results (xi, u, v) at the last time step to a json file (to be used as ground truth for inversion)
    xi_truth = vec(forward_simulation_results)[end][1:my_mesh_2D.numOfCells]
    wse_truth = xi_truth .+ wstill
    h_truth = xi_truth .+ hstill
    q_x_truth = vec(forward_simulation_results)[end][my_mesh_2D.numOfCells+1:2*my_mesh_2D.numOfCells]
    q_y_truth = vec(forward_simulation_results)[end][2*my_mesh_2D.numOfCells+1:3*my_mesh_2D.numOfCells]
    u_truth = q_x_truth ./ (h_truth .+ swe2d_extra_params.swe_2D_constants.h_small)
    v_truth = q_y_truth ./ (h_truth .+ swe2d_extra_params.swe_2D_constants.h_small)

    #If ManningN_option is variable_as_function_of_h, update ManningN_cell based on the ManningN_function_type and ManningN_function_parameters
    #if settings.forward_settings.ManningN_option == "variable_as_function_of_h"
    #    ManningN_cells = update_ManningN_forward_simulation(h_truth, settings)
    #end

    #compute bed slope 
    zb_ghostCells, zb_faces, S0_cells, S0_faces = update_bed_data(my_mesh_2D, zb_cells)

    #compute the friction terms
    friction_x_truth, friction_y_truth = compute_friction_terms(settings, h_truth, q_x_truth, q_y_truth, swe2d_extra_params.ManningN_cells, 
                 nothing, swe2d_extra_params, my_mesh_2D, swe2d_extra_params.swe_2D_constants.g, swe2d_extra_params.swe_2D_constants.k_n,
                  swe2d_extra_params.swe_2D_constants.h_small)

    # Function to replace NaN with nothing in arrays
    #replace_nan(x) = map(v -> isnan(v) ? nothing : v, x)
    replace_nan(x) = x

    #save the simulation results as truth to a json file
    open(joinpath(case_path, settings.forward_settings.save_solution_truth_file_name), "w") do io
        JSON3.pretty(io, Dict(
            "xi_truth" => replace_nan(xi_truth),
            "wstill_truth" => replace_nan(wstill),
            "hstill_truth" => replace_nan(hstill),
            "wse_truth" => replace_nan(wse_truth),
            "h_truth" => replace_nan(h_truth),
            "u_truth" => replace_nan(u_truth),
            "v_truth" => replace_nan(v_truth),
            "zb_cell_truth" => replace_nan(zb_cell_truth),
            "S0_cells_truth" => replace_nan(S0_cells),
            "S0_faces_truth" => replace_nan(S0_faces),
            "ManningN_cells_truth" => replace_nan(ManningN_cells),
            "friction_x_truth" => replace_nan(friction_x_truth),
            "friction_y_truth" => replace_nan(friction_y_truth),
            "ManningN_zone_values_truth" => replace_nan(ManningN_zone_values_truth),
            "inlet_discharges_truth" => replace_nan(inlet_discharges_truth)))
        println(io)
    end

    #save the simulation results to vtk files
    swe_2D_save_results_SciML(swe2d_extra_params, forward_simulation_results, friction_x_truth, friction_y_truth, S0_cells, S0_faces)
    
end


