#Some tools for 2D shallow water equations solver

#update scalar at ghost cells: ghost cells have the same value as the neighbor cell
function update_ghost_cells_scalar!(numOfAllBounaryFaces, allBoundaryFacesIDs_List, faceCells_Dict, scalar_cells, scalar_ghostCells)
    for iBoundaryFace in 1:numOfAllBounaryFaces
    
        if length(faceCells_Dict[allBoundaryFacesIDs_List[iBoundaryFace]]) !=1
            println("Error: the number of cells for boundary face ", allBoundaryFacesIDs_List[iBoundaryFace], " is not 1.")
            readline()
            exit(0)
        end
    
        cellID_neighbor = faceCells_Dict[allBoundaryFacesIDs_List[iBoundaryFace]][1]
        scalar_ghostCells[iBoundaryFace] = scalar_cells[cellID_neighbor]
    end
    #println("scalar_ghostCells: ", scalar_ghostCells)
end

# computer gradient of a scalar field
function compute_scalar_gradients!(numOfCells, cell_areas, cell_normals, face_lengths, cellNodesCount, cellFacesList, cellNeighbors_Dict, scalar_variable, grad_scalar_variable)
    #check the size of the grad_scalar_variable
    if size(grad_scalar_variable) != (numOfCells, 2)
        println("Error: grad_scalar_variable size is not correct.")
        readline()
        exit(-1)
    end

    fill!(grad_scalar_variable, 0.0)
    
    #loop over all cells to compute the gradient of the scalar field
    for iCell in 1:numOfCells
        cell_gradient = [0.0, 0.0]  # Gradient accumulator for this cell
        
        #neighbor cells of the current cell
        cellNeighbors = cellNeighbors_Dict[iCell]
        
        #number of nodes for the current cell
        nNodes = cellNodesCount[iCell]
        
        cell_faces = cellFacesList[iCell,:]
        
        #loop over all faces of the current cell
        for iFace in 1:nNodes
            faceID = cell_faces[iFace]
            neighbor_cellID = cellNeighbors[iFace]
            
            # Value of the variable at the current cell and neighbor cell
            variable_c = scalar_variable[iCell]
            if neighbor_cellID < 0  # Boundary face
                variable_n = variable_c  # Assume zero gradient at boundary
            else
                variable_n = scalar_variable[neighbor_cellID]
            end
            
            # Compute the variable value on face (average of cell and neighbor for now; can be improved with interpolation)
            variable_f = (variable_c + variable_n) / 2.0
            
            # Compute flux contribution
            flux = cell_normals[iCell][iFace] * variable_f * face_lengths[abs(faceID)]
            cell_gradient .+= flux
        end
        
        # Finalize the gradient by dividing by the cell area
        grad_scalar_variable[iCell, :] = cell_gradient / cell_areas[iCell]
        
    end
end

# interpolate a scalar field from cell centers to face centers
function cells_to_faces_scalar!(numOfFaces, faceCells_Dict, scalar_variable_c, scalar_variable_f)
    
    #loop through faces  
    @inbounds for iFace in 1:numOfFaces
        #get current face's cells
        facecells = faceCells_Dict[iFace]
        
        if length(facecells) == 2   #internal face
            scalar_variable_f[iFace] = (scalar_variable_c[facecells[1]] + scalar_variable_c[facecells[2]]) / 2.0    
        else  #boundary face
            scalar_variable_f[iFace] = scalar_variable_c[facecells[1]]  
        end
    end
end



