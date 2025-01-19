
#update scalar at ghost cells: ghost cells have the same value as the neighbor cell
function update_ghost_cells_scalar(my_mesh_2D::mesh_2D, scalar_cells::AbstractVector{T})::AbstractVector{T} where T<:Real

    # Check if scalar_cells is an AbstractArray
    #@assert isa(scalar_cells, AbstractArray) "scalar_cells must be an AbstractArray"

    # Check if scalar_cells can handle dual numbers
    #@assert eltype(scalar_cells) <: Real "scalar_cells must have elements of a real number type"

    scalar_ghostCells = [
        let faceID = my_mesh_2D.allBoundaryFacesIDs_List[iBoundaryFace]
            faceCells = my_mesh_2D.faceCells_Dict[faceID]
            
            # Check validity
            #length(faceCells) == 1 || throw(error("Error: the number of cells for boundary face $(faceID) is not 1."))

            #@assert length(faceCells) == 1 "Error: the number of cells for boundary face $(faceID) is not 1."

            # Ensure index is within bounds
            #@assert 1 <= faceCells[1] <= length(scalar_cells) "Index out of bounds for scalar_cells"
            
            # Get value from neighboring cell (for boundary faces, the neighboring cell is the same as the current cell)
            scalar_cells[faceCells[1]]
        end
        for iBoundaryFace in 1:my_mesh_2D.numOfAllBounaryFaces
    ]

    return scalar_ghostCells
end

# computer gradient of a scalar field
function compute_scalar_gradients_from_cells(my_mesh_2D::mesh_2D, scalar_variable::AbstractVector{T})::Matrix{T} where T<:Real

    # Construct the gradient matrix directly
    grad_scalar_variable = [compute_cell_gradient_from_cells(iCell, my_mesh_2D, scalar_variable) 
                            for iCell in 1:my_mesh_2D.numOfCells]

    gradients = [grad_scalar_variable[i][j] for i in eachindex(grad_scalar_variable), j in eachindex(grad_scalar_variable[1])]
     
    return gradients
end

function compute_cell_gradient_from_cells(iCell::Int, my_mesh_2D::mesh_2D, scalar_variable::AbstractVector{T})::Vector{T} where T<:Real
    # ... cell gradient computation logic ...
    cell_gradient = [zero(T), zero(T)]  # Gradient accumulator for this cell
    #cell_gradient = SVector{2,T}(zero(T), zero(T))
        
    #neighbor cells of the current cell
    cellNeighbors = my_mesh_2D.cellNeighbors_Dict[iCell]
    
    #number of nodes for the current cell
    nNodes = my_mesh_2D.cellNodesCount[iCell]
    
    cell_faces = my_mesh_2D.cellFacesList[iCell,:]

    #@show iCell
    #@show cellNeighbors
    #@show cell_faces

    #loop over all faces of the current cell
    for iFace in 1:nNodes
        faceID = cell_faces[iFace]
        neighbor_cellID = cellNeighbors[iFace]
        
        # Value of the variable at the current cell and neighbor cell
        variable_c = scalar_variable[iCell]
        if my_mesh_2D.bFace_is_boundary[faceID]  # If it is a boundary face
            variable_n = variable_c  # Assume zero gradient at boundary
        else
            variable_n = scalar_variable[neighbor_cellID]
        end

        # Compute the variable value on face (average of cell and neighbor for now; can be improved with interpolation)
        variable_f = T((variable_c + variable_n) / 2.0)
        
        # Compute flux contribution
        flux_temp = my_mesh_2D.cell_normals[iCell][iFace] * variable_f * my_mesh_2D.face_lengths[faceID]
        
        cell_gradient = cell_gradient + flux_temp

        #@show iCell, iFace, faceID, neighbor_cellID
        #@show variable_c
        #@show variable_n
        #@show my_mesh_2D.cell_normals[iCell][iFace]
        #@show flux_temp
        #@show cell_gradient

    end
    
    return cell_gradient / my_mesh_2D.cell_areas[iCell]
end

# interpolate a scalar field from cell centers to face centers
# The boundary condition is not considered here.
function cells_to_faces_scalar(my_mesh_2D::mesh_2D, scalar_variable_c::AbstractVector{T})::Vector{T} where T<:Real

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

#interpolate a scalar field from nodes to cells: cell values are the average of the node values
function nodes_to_cells_scalar(my_mesh_2D::mesh_2D, scalar_variable_n::AbstractVector{T})::Vector{T} where T<:Real

    scalar_variable_c = [
        mean(scalar_variable_n[my_mesh_2D.cellNodesList[iCell,:][1:my_mesh_2D.cellNodesCount[iCell]]])
        for iCell in 1:my_mesh_2D.numOfCells
    ]

    return scalar_variable_c
end

# interploate a scalar field from nodes to face centers
# Note: faceNodes_r_Dict is the reverse of faceNodes_Dict. faceNodes_r_Dict[iFace] = (node1, node2)
function nodes_to_faces_scalar(my_mesh_2D::mesh_2D, scalar_variable_n::AbstractVector{T})::Vector{T} where T<:Real

    scalar_variable_f = [
        (scalar_variable_n[my_mesh_2D.faceNodes_r_Dict[iFace][1]] + 
         scalar_variable_n[my_mesh_2D.faceNodes_r_Dict[iFace][2]]) / 2.0
        for iFace in 1:my_mesh_2D.numOfFaces
    ]

    return scalar_variable_f
end

# computer gradient of a scalar field from values at faces
# Before calling this function, the scalar field should be interpolated from cells to faces.
function compute_scalar_gradients_from_faces(my_mesh_2D::mesh_2D, scalar_variable_f::AbstractVector{T})::Matrix{T} where T<:Real

    # Construct the gradient matrix directly
    grad_scalar_variable = [compute_cell_gradient_from_faces(iCell, my_mesh_2D, scalar_variable_f) 
                            for iCell in 1:my_mesh_2D.numOfCells]

    gradients = [grad_scalar_variable[i][j] for i in eachindex(grad_scalar_variable), j in eachindex(grad_scalar_variable[1])]
   
         
    return gradients
end

function compute_cell_gradient_from_faces(iCell::Int, my_mesh_2D::mesh_2D, scalar_variable_f::AbstractVector{T})::Vector{T} where T<:Real
    # ... cell gradient computation logic ...
    cell_gradient = [zero(T), zero(T)]  # Gradient accumulator for this cell    
    #cell_gradient = @SVector [zero(T), zero(T)]
   
    #number of nodes for the current cell
    nNodes = my_mesh_2D.cellNodesCount[iCell]
    
    cell_faces = my_mesh_2D.cellFacesList[iCell,:]

    #loop over all faces of the current cell
    for iFace in 1:nNodes
        faceID = cell_faces[iFace]

        # Compute flux contribution
        flux_temp = my_mesh_2D.cell_normals[iCell][iFace] * scalar_variable_f[faceID] * my_mesh_2D.face_lengths[faceID]
        
        cell_gradient = cell_gradient + flux_temp

        if iCell == 1065
            @show iCell, iFace, faceID
            @show my_mesh_2D.cell_normals[iCell][iFace]
            @show scalar_variable_f[faceID]
            @show my_mesh_2D.face_lengths[faceID]
            @show flux_temp
            @show cell_gradient
        end

        #@show iCell, iFace, faceID
        #@show my_mesh_2D.cell_normals[iCell][iFace]
        #@show flux_temp
        #@show cell_gradient

    end

    if iCell == 1065
        @show cell_gradient / my_mesh_2D.cell_areas[iCell]
    end
    
    return cell_gradient / my_mesh_2D.cell_areas[iCell]
end




