
#update scalar at ghost cells: ghost cells have the same value as the neighbor cell
function update_ghost_cells_scalar(my_mesh_2D, scalar_cells)

    # Check if scalar_cells is an AbstractArray
    @assert isa(scalar_cells, AbstractArray) "scalar_cells must be an AbstractArray"

    # Check if scalar_cells can handle dual numbers
    @assert eltype(scalar_cells) <: Real "scalar_cells must have elements of a real number type"

    scalar_ghostCells = [
        let faceID = my_mesh_2D.allBoundaryFacesIDs_List[iBoundaryFace]
            faceCells = my_mesh_2D.faceCells_Dict[faceID]
            
            # Check validity
            #length(faceCells) == 1 || throw(error("Error: the number of cells for boundary face $(faceID) is not 1."))

            @assert length(faceCells) == 1 "Error: the number of cells for boundary face $(faceID) is not 1."

            # Ensure index is within bounds
            @assert 1 <= faceCells[1] <= length(scalar_cells) "Index out of bounds for scalar_cells"
            
            # Get value from neighboring cell (for boundary faces, the neighboring cell is the same as the current cell)
            scalar_cells[faceCells[1]]
        end
        for iBoundaryFace in 1:my_mesh_2D.numOfAllBounaryFaces
    ]

    return scalar_ghostCells
end

# computer gradient of a scalar field
function compute_scalar_gradients(my_mesh_2D, scalar_variable)

    # Construct the gradient matrix directly
    grad_scalar_variable = [compute_cell_gradient(iCell, my_mesh_2D, scalar_variable) 
                            for iCell in 1:my_mesh_2D.numOfCells]

    gradients = [grad_scalar_variable[i][j] for i in eachindex(grad_scalar_variable), j in eachindex(grad_scalar_variable[1])]
     
    return gradients
end

function compute_cell_gradient(iCell, my_mesh_2D, scalar_variable)
    # ... cell gradient computation logic ...
    cell_gradient = eltype(scalar_variable).([0.0, 0.0])  # Gradient accumulator for this cell
        
    #neighbor cells of the current cell
    cellNeighbors = my_mesh_2D.cellNeighbors_Dict[iCell]
    
    #number of nodes for the current cell
    nNodes = my_mesh_2D.cellNodesCount[iCell]
    
    cell_faces = my_mesh_2D.cellFacesList[iCell,:]
    
    #loop over all faces of the current cell
    for iFace in 1:nNodes
        faceID = cell_faces[iFace]
        neighbor_cellID = cellNeighbors[iFace]
        
        # Value of the variable at the current cell and neighbor cell
        variable_c = scalar_variable[iCell]
        if neighbor_cellID > my_mesh_2D.numOfCells  # Boundary face if the neighbor cell ID is greater than the number of cells
            variable_n = variable_c  # Assume zero gradient at boundary
        else
            variable_n = scalar_variable[neighbor_cellID]
        end
        
        # Compute the variable value on face (average of cell and neighbor for now; can be improved with interpolation)
        variable_f = (variable_c + variable_n) / 2.0
        
        # Compute flux contribution
        flux_temp = my_mesh_2D.cell_normals[iCell][iFace] * variable_f * my_mesh_2D.face_lengths[faceID]
        
        cell_gradient = cell_gradient + flux_temp

    end
    
    return cell_gradient / my_mesh_2D.cell_areas[iCell]
end

# interpolate a scalar field from cell centers to face centers
# The boundary condition is not considered here.
function cells_to_faces_scalar(my_mesh_2D, scalar_variable_c)
    # # Create new array instead of modifying in-place
    # scalar_variable_f = zeros(eltype(scalar_variable_c), my_mesh_2D.numOfFaces)
    
    # #loop through faces  
    # for iFace in 1:my_mesh_2D.numOfFaces
    #     #get current face's cells
    #     facecells = my_mesh_2D.faceCells_Dict[iFace]
        
    #     if length(facecells) == 2   #internal face
    #         scalar_variable_f[iFace] = (scalar_variable_c[facecells[1]] + scalar_variable_c[facecells[2]]) / 2.0    
    #     else  #boundary face
    #         scalar_variable_f[iFace] = scalar_variable_c[facecells[1]]  
    #     end
    # end

    scalar_variable_f = [
        if length(my_mesh_2D.faceCells_Dict[iFace]) == 2
            # Internal face - average of two cells
            (scalar_variable_c[my_mesh_2D.faceCells_Dict[iFace][1]] + 
             scalar_variable_c[my_mesh_2D.faceCells_Dict[iFace][2]]) / 2.0
        else
            # Boundary face - use the single cell value
            scalar_variable_c[my_mesh_2D.faceCells_Dict[iFace][1]]
        end
        for iFace in 1:my_mesh_2D.numOfFaces
    ]

    return scalar_variable_f
    
end



