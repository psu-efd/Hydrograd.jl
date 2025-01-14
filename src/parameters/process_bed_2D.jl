# process bed profiles, such as setup bed profile, interplation of bed elevation from face to cell,
# cell to face, calculate slope S0, etc. 
#
# This file should be problem specific (each problem may have different bed profile)

#Set up the bed elevation: assuming the nodes coordinates already have the bed elevation
#   interpolate the bed elevation from nodes to cells, then cells to faces, then compute the bed slope at cells
#   
function setup_bed(settings::ControlSettings, my_mesh_2D::mesh_2D, nodeCoordinates::Matrix{T}, case_path::String, bPlotBed::Bool=false) where {T<:Real}
    zb_nodes = @view nodeCoordinates[:, 3]

    # Initialize bed elevation at cell centers
    #zb_cells = zeros(Float64, my_mesh_2D.numOfCells)  # zb at cell centers 

    #update zb at cells, ghost cells, faces, and S0 from nodeCoordinates
    #zb_cells, zb_ghostCells, zb_faces, S0 = interploate_zb_from_nodes_to_cells_ghostCells_faces_and_compute_S0(my_mesh_2D, nodeCoordinates[:, 3])

    #interpolate zb from nodes to cells
    zb_cells = nodes_to_cells_scalar(my_mesh_2D, zb_nodes)
   
    #update the ghost cells for zb from zb_cells
    zb_ghostCells = update_ghost_cells_scalar(my_mesh_2D, zb_cells)

    #interpolate zb from nodes to faces
    #zb_faces = nodes_to_faces_scalar(my_mesh_2D, zb_nodes)
    #or interpolate zb from cells to faces
    zb_faces = cells_to_faces_scalar(my_mesh_2D, zb_cells)

    #compute the bed slope at cells from zb_faces
    S0 = -1.0 * compute_scalar_gradients_from_faces(my_mesh_2D, zb_faces)

    #@show typeof(S0)

    # optionally plot the bed for checking
    if bPlotBed
        vector_data = [S0]
        vector_names = ["S0"]

        scalar_data = [zb_cells]
        scalar_names = ["zb_cells"]

        file_path = joinpath(case_path, "zb_S0.vtk")
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount,
            "", "float", 0.0, scalar_data, scalar_names, vector_data, vector_names)

        if settings.bVerbose
            println("zb and S0 are saved to ", file_path)
        end
    end

    return zb_cells, zb_ghostCells, zb_faces, S0

end

function interploate_zb_from_cells_to_ghostCells_faces_and_compute_S0(my_mesh_2D::mesh_2D, zb_cells::AbstractVector{T}) where {T<:Real}
   
    #update the ghost cells for zb from zb_cells
    zb_ghostCells = update_ghost_cells_scalar(my_mesh_2D, zb_cells)

    #interpolate zb from cells to faces
    zb_faces = cells_to_faces_scalar(my_mesh_2D, zb_cells)

    #compute the bed slope at cells from zb_faces
    S0 = -1.0 * compute_scalar_gradients_from_faces(my_mesh_2D, zb_faces)

    #@show typeof(S0)
    #@show size(S0)
    #@show S0[1]
    
    return zb_ghostCells, zb_faces, S0

end

#interpolate zb from nodes to face and compute the bed slope at cells
# function interploate_zb_from_nodes_to_cells_ghostCells_faces_and_compute_S0(my_mesh_2D::mesh_2D, zb_nodes::AbstractVector{T}) where {T<:Real}

#     #interpolate zb from nodes to cells
#     zb_cells = nodes_to_cells_scalar(my_mesh_2D, zb_nodes)
   
#     #update the ghost cells for zb from zb_cells
#     zb_ghostCells = update_ghost_cells_scalar(my_mesh_2D, zb_cells)

#     #interpolate zb from nodes to faces
#     #zb_faces = nodes_to_faces_scalar(my_mesh_2D, zb_nodes)

#     #compute the bed slope at cells from zb_faces
#     S0 = -1.0 * compute_scalar_gradients_from_faces(my_mesh_2D, zb_faces)
    
#     return zb_cells, zb_ghostCells, zb_faces, S0

# end
