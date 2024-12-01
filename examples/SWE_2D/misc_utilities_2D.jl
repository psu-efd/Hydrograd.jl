#Some tools for 2D shallow water equations solver

#update scalar at ghost cells: ghost cells have the same value as the neighbor cell
function update_ghost_cells_scalar!(numOfAllBounaryFaces, allBoundaryFacesIDs_List, faceCells_Dict, scalar_cells, scalar_ghostCells)
    for iBoundaryFace in 1:numOfAllBounaryFaces
    
        if length(faceCells_Dict[allBoundaryFacesIDs_List[iBoundaryFace]]) !=1
            println("Error: the number of cells for boundary face ", allBoundaryFacesIDs_List[iBoundaryFace], " is not 1.")
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


function export_to_vtk(filename, nodeCoordinates, cellNodesList, cellNodesCount, scalar_data, scalar_names, vector_data, vector_names)
    # Open the file for writing
    open(filename, "w") do file
        # Write the VTK header
        println(file, "# vtk DataFile Version 2.0")
        println(file, "2D Unstructured Mesh")
        println(file, "ASCII")
        println(file, "DATASET UNSTRUCTURED_GRID")
        
        # Write the nodes
        num_nodes = size(nodeCoordinates, 1)
        println(file, "POINTS $num_nodes float")
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
            println(file, "SCALARS $name float 1")
            println(file, "LOOKUP_TABLE default")
            for value in scalar
                println(file, value)
            end
        end
        
        # Write vector data
        for (vector, name) in zip(vector_data, vector_names)
            println(file, "VECTORS $name float")
            for vec in eachrow(vector)
                println(file, "$(vec[1]) $(vec[2]) 0.0")  # Add 0.0 for the z-component
            end
        end
    end
end

