#Some tools for 2D shallow water equations solver

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