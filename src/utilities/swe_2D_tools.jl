#some tools for 2D SWE

#calculate the total water volume in the domain
function swe_2D_calc_total_water_volume(h, my_mesh_2D)
    total_water_volume = sum(h .* my_mesh_2D.cell_areas)
    return total_water_volume
end

#save results (sol is a solution from SciML ODE solver)
function swe_2D_save_results_SciML(swe2d_extra_params, sol, friction_x_truth, friction_y_truth, S0_cells, S0_faces)
    
    settings = swe2d_extra_params.settings
    my_mesh_2D = swe2d_extra_params.my_mesh_2D
    nodeCoordinates = swe2d_extra_params.nodeCoordinates
    zb_cells = swe2d_extra_params.zb_cells
    save_path = swe2d_extra_params.case_path

    ManningN_cells = swe2d_extra_params.ManningN_cells
    ks_cells = swe2d_extra_params.ks_cells
    h_ks_cells = swe2d_extra_params.h_ks_cells
    friction_factor_cells = swe2d_extra_params.friction_factor_cells
    Re_cells = swe2d_extra_params.Re_cells

    wstill = swe2d_extra_params.wstill    
    hstill = swe2d_extra_params.hstill

    total_water_volume = []

    #calculate total volume of water 
    for (index, state) in enumerate(sol)

        field_name = "forward_simulation_saved_index"
        field_type = "integer"
        field_value = index

        # Extract solution components
        xi_array = state[1:my_mesh_2D.numOfCells]
        h_array = xi_array .+ hstill
        q_x_array = state[my_mesh_2D.numOfCells+1:2*my_mesh_2D.numOfCells]
        q_y_array = state[2*my_mesh_2D.numOfCells+1:3*my_mesh_2D.numOfCells]
        u_array = q_x_array ./ (h_array .+ swe2d_extra_params.swe_2D_constants.h_small)
        v_array = q_y_array ./ (h_array .+ swe2d_extra_params.swe_2D_constants.h_small)
        Umag_array = sqrt.(u_array.^2 .+ v_array.^2)

        #If ManningN_option is variable, update ManningN_cell based on the ManningN_function_type and ManningN_function_parameters
        if settings.forward_settings.ManningN_option == "variable"
            ManningN_cells, h_ks_cells, friction_factor_cells, Re_cells = update_ManningN_forward_simulation(h_array, Umag_array, ks_cells, settings)
        end

        #computer wet/dry flags
        b_dry_wet, b_Adjacent_to_dry_land, b_Adjacent_to_high_dry_land = process_dry_wet_flags(my_mesh_2D, h_array, zb_cells, swe2d_extra_params.swe_2D_constants)

        #convert b_dry_wet, b_Adjacent_to_dry_land, b_Adjacent_to_high_dry_land to floats
        b_dry_wet_float = Float64.(b_dry_wet)
        b_Adjacent_to_dry_land_float = Float64.(b_Adjacent_to_dry_land)
        b_Adjacent_to_high_dry_land_float = Float64.(b_Adjacent_to_high_dry_land)

        #calculate total volume of water
        push!(total_water_volume, swe_2D_calc_total_water_volume(h_array, my_mesh_2D))

        u_temp = q_x_array ./ (h_array .+ swe2d_extra_params.swe_2D_constants.h_small)
        v_temp = q_y_array ./ (h_array .+ swe2d_extra_params.swe_2D_constants.h_small)
        U_vector = hcat(u_temp, v_temp)

        #add slope to the vector data
        slope_x = S0_cells[:,1]
        slope_y = S0_cells[:,2]
        slope_vector = hcat(slope_x, slope_y)
        
                        
        vector_data = [U_vector, slope_vector] 
        vector_names = ["U", "slope"]

        WSE = h_array + zb_cells
            
        scalar_data = [xi_array, wstill, hstill, h_array, q_x_array, q_y_array, ManningN_cells, ks_cells, h_ks_cells, friction_factor_cells, Re_cells, zb_cells, WSE, friction_x_truth, friction_y_truth, b_dry_wet_float, b_Adjacent_to_dry_land_float, b_Adjacent_to_high_dry_land_float]
        scalar_names = ["xi", "wstill", "hstill", "h", "hu", "hv", "ManningN", "ks", "h_ks", "friction_factor", "Re", "zb_cell", "WSE", "friction_x", "friction_y", "b_dry_wet", "b_Adjacent_to_dry_land", "b_Adjacent_to_high_dry_land"]

        vtk_fileName = @sprintf("forward_simulation_results_%04d.vtk", index)
            
        file_path = joinpath(save_path, vtk_fileName ) 
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, field_name, field_type, field_value, scalar_data, scalar_names, vector_data, vector_names)    
    
    end
        
    open(joinpath(save_path, "total_water_volume.csv"), "w") do fo
        println(fo, "total_water_volume")
        for volume in total_water_volume
            println(fo, volume)
        end
    end
end

#save results (sol is a solution from custom ODE solver)
function swe_2D_save_results_custom(sol, swe2d_extra_params)

    my_mesh_2D = swe2d_extra_params.my_mesh_2D
    nodeCoordinates = swe2d_extra_params.nodeCoordinates
    zb_cells = swe2d_extra_params.zb_cells
    save_path = swe2d_extra_params.case_path

    #calculate total volume of water 
    total_water_volume = []

    #for index in 1:size(sol, 3)
    for index in axes(sol, 2)

        field_name = "forward_simulation_saved_index"
        field_type = "integer"
        field_value = index

        Q = @view sol[:,index]  # Q will be a view of the 2D slice at timestep index

        push!(total_water_volume, swe_2D_calc_total_water_volume(Q[1:my_mesh_2D.numOfCells], my_mesh_2D))

        u_temp = Q[my_mesh_2D.numOfCells+1:2*my_mesh_2D.numOfCells] ./ Q[1:my_mesh_2D.numOfCells]
        v_temp = Q[2*my_mesh_2D.numOfCells+1:3*my_mesh_2D.numOfCells] ./ Q[1:my_mesh_2D.numOfCells]
        U_vector = hcat(u_temp, v_temp)
                        
        vector_data = [U_vector] 
        vector_names = ["U"]

        WSE = Q[1:my_mesh_2D.numOfCells] + zb_cells
            
        scalar_data = [Q[1:my_mesh_2D.numOfCells], Q[my_mesh_2D.numOfCells+1:2*my_mesh_2D.numOfCells], Q[2*my_mesh_2D.numOfCells+1:3*my_mesh_2D.numOfCells], zb_cells, WSE]
        scalar_names = ["h", "hu", "hv", "zb_cell", "WSE"]

        vtk_fileName = @sprintf("forward_simulation_results_%04d.vtk", index)
            
        file_path = joinpath(save_path, vtk_fileName ) 
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, field_name, field_type, field_value, scalar_data, scalar_names, vector_data, vector_names)    
    
    end

        
    # open(joinpath(save_path, "total_water_volume.csv"), "w") do fo
    #     println(fo, "total_water_volume")
    #     for volume in total_water_volume
    #         println(fo, volume)
    #     end
    # end
end


function export_to_vtk_2D(filename, nodeCoordinates, cellNodesList, cellNodesCount, field_name, field_type, field_value, scalar_data, scalar_names, vector_data, vector_names)

    #check if field_name, field_type, field_value are valid
    if !isa(field_name, String) || !isa(field_type, String) || !isa(field_value, Number)
        println("field_name, field_type, field_value are not valid")
        println("field_name: ", field_name)
        println("field_type: ", field_type)
        println("field_value: ", field_value)

        return
    end

    # Open the file for writing
    open(filename, "w") do file
        # Write the VTK header
        println(file, "# vtk DataFile Version 2.0")
        println(file, "2D Unstructured Mesh")
        println(file, "ASCII")
        println(file, "DATASET UNSTRUCTURED_GRID")

        # Write the field name and value (if it is not empty)
        if field_name != ""
            println(file, "FIELD FieldData 1")
            println(file, "$field_name 1 1 $field_type")
            println(file, "$field_value")
        end
        
        # Write the nodes
        num_nodes = size(nodeCoordinates, 1)
        println(file, "POINTS $num_nodes double")
        for row in eachrow(nodeCoordinates)
            println(file, "$(row[1]) $(row[2]) $(row[3])")
        end
        
        # Write the cells
        num_cells = size(cellNodesList, 1)
        total_cell_size = sum(nNodes + 1 for nNodes in cellNodesCount)  # Include the cell type
        println(file, "CELLS $num_cells $total_cell_size")
        for cellID in 1:num_cells
            cell = cellNodesList[cellID, :][1:cellNodesCount[cellID]]
            println(file, "$(length(cell)) $(join(cell .- 1, ' '))")  # VTK uses 0-based indexing
        end
        
        # Write the cell types (assume all are polygons for 2D unstructured meshes)
        println(file, "CELL_TYPES $num_cells")
        for _ in 1:num_cells
            println(file, "7")  # 7 corresponds to VTK_POLYGON
        end
        
        # Write cell data
        println(file, "CELL_DATA $num_cells")
        
        # Write scalar data
        for (scalar, name) in zip(scalar_data, scalar_names)
            println(file, "SCALARS $name double 1")
            println(file, "LOOKUP_TABLE default")
            for value in scalar
                println(file, value)
            end
        end
        
        # Write vector data
        for (vector, name) in zip(vector_data, vector_names)
            println(file, "VECTORS $name double")
            for vec in eachrow(vector)
                println(file, "$(vec[1]) $(vec[2]) 0.0")  # Add 0.0 for the z-component
            end
        end
    end
end



