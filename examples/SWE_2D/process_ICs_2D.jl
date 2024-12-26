# process initial condition, such as setting up initial free surface, water depth, etc. 
# This file should be problem specific because each problem should have different ICs. 

# setup initial condition: wse, h, q_x, q_y
function setup_initial_condition!(forward_simulation_initial_condition_options, forward_simulation_initial_condition_constant_values, 
    forward_simulation_initial_condition_values_from_file, my_mesh_2D, nodeCoordinates, eta, zb_cell, h, q_x, q_y,swe_2D_constants, bPlot::Bool=false)
    
    if forward_simulation_initial_condition_options == "constant"
        #loop over cells
        @inbounds for i in 1:my_mesh_2D.numOfCells
            h[i] = forward_simulation_initial_condition_constant_values[1]
            q_x[i] = forward_simulation_initial_condition_constant_values[2]
            q_y[i] = forward_simulation_initial_condition_constant_values[3]
        end
    elseif forward_simulation_initial_condition_options == "from_file"
        #make sure the length of the initial condition values is the same as the number of cells
        if length(forward_simulation_initial_condition_values_from_file["h"]) != my_mesh_2D.numOfCells
            error("The length of the initial condition values of h is not the same as the number of cells.")
        end
        if length(forward_simulation_initial_condition_values_from_file["q_x"]) != my_mesh_2D.numOfCells
            error("The length of the initial condition values of q_x is not the same as the number of cells.")
        end
        if length(forward_simulation_initial_condition_values_from_file["q_y"]) != my_mesh_2D.numOfCells
            error("The length of the initial condition values of q_y is not the same as the number of cells.")
        end

        #copy the initial condition values
        h = deepcopy(forward_simulation_initial_condition_values_from_file["h"])
        q_x = deepcopy(forward_simulation_initial_condition_values_from_file["q_x"])
        q_y = deepcopy(forward_simulation_initial_condition_values_from_file["q_y"])
    end
    
    h[h.<0.0] .= swe_2D_constants.h_small  #ensure positivity of water depth h
    
    #update the free surface elevation again in case h has been clipped
    @inbounds for i in 1:my_mesh_2D.numOfCells
        eta[i] = h[i] + zb_cell[i]
    end
  
    #optionally plot the ICs for checking 
    if bPlot
        vector_data = [] 
        vector_names = []
        
        scalar_data = [eta, h, q_x, q_y, zb_cell]
        scalar_names = ["eta", "h", "q_x", "q_y", "zb_cell"]
        
        file_path = joinpath(@__DIR__, "initial_conditions.vtk" ) 
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, 
                         scalar_data, scalar_names, vector_data, vector_names)    
        println("initial conditions are saved to ", file_path)
    end
    
end

#update ghost cells for eta, h, q_x, q_y
function setup_ghost_cells_initial_condition!(my_mesh_2D, eta, h, q_x, q_y, eta_ghostCells, h_ghostCells, q_x_ghostCells, q_y_ghostCells)

    for iBoundaryFace in 1:my_mesh_2D.numOfAllBounaryFaces
        cellID_neighbor = my_mesh_2D.faceCells_Dict[my_mesh_2D.allBoundaryFacesIDs_List[iBoundaryFace]][1]
        
        eta_ghostCells[iBoundaryFace] = eta[cellID_neighbor]
        h_ghostCells[iBoundaryFace] = h[cellID_neighbor]
        q_x_ghostCells[iBoundaryFace] = q_x[cellID_neighbor]
        q_y_ghostCells[iBoundaryFace] = q_y[cellID_neighbor]
    end
    
    #println("h_ghostCells: ", h_ghostCells)
end 


