# process bed profiles, such as setup bed profile, interplation of bed elevation from face to cell,
# cell to face, calculate slope S0, etc. 
#
# This file should be problem specific (each problem may have different bed profile)



#loss due to slope regularization
function calc_slope_loss(zb_cells, my_mesh_2D)
    
     # Compute bed slope at cell centers (Bed slope is the negative of zb gradient)
     S0 = -1.0 * compute_scalar_gradients(my_mesh_2D, zb_cells)

     #compute the magnitude of S0
     S0_mag = sqrt.(S0[:,1].^2 + S0[:,2].^2)  

    S0_max = 0.2
    S0_center = 0.0

    #return max.(S0_center, (slope_temp .- S0_max)) 
    return sum(max.(S0_center, (S0_mag .- S0_max)))
end

#Set up the bed elevation: assuming the nodes coordinates already have the bed elevation
#   interpolate the bed elevation from nodes to cells, then cells to faces, then compute the bed slope at cells
#   
function setup_bed(settings, my_mesh_2D, nodeCoordinates, bPlotBed::Bool=false)

     # Initialize bed elevation at cell centers
     zb_cells = zeros(Float64, my_mesh_2D.numOfCells)  # zb at cell centers 

    #loop through cells
    for i in 1:my_mesh_2D.numOfCells
        #get the node coordinates of the cell
        cell_nodes = my_mesh_2D.cellNodesList[i,:][1:my_mesh_2D.cellNodesCount[i]]
        cell_node_coordinates = nodeCoordinates[cell_nodes, :]

        #interpolate the bed elevation from nodes to cell centroid
        zb_cells[i] = mean(cell_node_coordinates[:, 3])
    end

    #interpolate zb from cell to face and compute the bed slope at cells
    zb_ghostCells, zb_faces, S0 = interploate_zb_from_cell_to_face_and_compute_S0(my_mesh_2D, zb_cells)

    # optionally plot the bed for checking
    if bPlotBed
        vector_data = [S0] 
        vector_names = ["S0"]

        scalar_data = [zb_cells]
        scalar_names = ["zb_cells"]
        
        file_path = joinpath(@__DIR__, "zb_S0.vtk" ) 
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, 
                         scalar_data, scalar_names, vector_data, vector_names)    
        
        if settings.bVerbose
            println("zb and S0 are saved to ", file_path)
        end
    end

    return zb_cells, zb_ghostCells, zb_faces, S0
    
end

#interpolate zb from cell to face and compute the bed slope at cells
function interploate_zb_from_cell_to_face_and_compute_S0(my_mesh_2D, zb_cells)

     # Setup zb_ghostCells
     zb_ghostCells = update_ghost_cells_scalar(my_mesh_2D, zb_cells)
     
     # Interpolate zb from cell to face
     zb_faces = cells_to_faces_scalar(my_mesh_2D, zb_cells)
     
     # Compute bed slope at cell centers (Bed slope is the negative of zb gradient)
     S0 = -1.0 * compute_scalar_gradients(my_mesh_2D, zb_cells)
          
     return zb_ghostCells, zb_faces, S0
    
end


