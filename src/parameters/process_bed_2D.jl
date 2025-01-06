# process bed profiles, such as setup bed profile, interplation of bed elevation from face to cell,
# cell to face, calculate slope S0, etc. 
#
# This file should be problem specific (each problem may have different bed profile)

#Set up the bed elevation: assuming the nodes coordinates already have the bed elevation
#   interpolate the bed elevation from nodes to cells, then cells to faces, then compute the bed slope at cells
#   
function setup_bed(settings::ControlSettings, my_mesh_2D::mesh_2D, nodeCoordinates::Matrix{T}, case_path::String, bPlotBed::Bool=false) where T <: Real

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
        
        file_path = joinpath(case_path, "zb_S0.vtk" ) 
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, 
                         "",  "float", 0.0, scalar_data, scalar_names, vector_data, vector_names)    
        
        if settings.bVerbose
            println("zb and S0 are saved to ", file_path)
        end
    end

    return zb_cells, zb_ghostCells, zb_faces, S0
    
end

#interpolate zb from cell to face and compute the bed slope at cells
function interploate_zb_from_cell_to_face_and_compute_S0(my_mesh_2D::mesh_2D, zb_cells::AbstractVector{T}) where T <: Real

     # Setup zb_ghostCells
     zb_ghostCells = update_ghost_cells_scalar(my_mesh_2D, zb_cells)
     
     # Interpolate zb from cell to face
     zb_faces = cells_to_faces_scalar(my_mesh_2D, zb_cells)
     
     # Compute bed slope at cell centers (Bed slope is the negative of zb gradient)
     S0 = -1.0 * compute_scalar_gradients(my_mesh_2D, zb_cells)

     #@show typeof(S0)
     #@show size(S0)
     #@show S0
          
     return zb_ghostCells, zb_faces, S0
    
end


