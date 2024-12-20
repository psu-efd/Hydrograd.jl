#some tools for 1D SWE

using Printf
using DifferentialEquations
using CSV 
using DataFrames
using Plots

#calculate the total water volume in the domain
function swe_2D_calc_total_water_volume(h, my_mesh_2D)
    total_water_volume = sum(h .* my_mesh_2D.cell_areas)
    return total_water_volume
end

#save results 
function swe_2D_save_results(sol, total_water_volume, my_mesh_2D, zb_cell, save_path)

    #solution at the end of the simulation
    #Q_final = sol.u[end]

    #calculate total volume of water 
    for (index, item) in enumerate(tuples(sol))

        Q = item[1]

        push!(total_water_volume, swe_2D_calc_total_water_volume(Q[:, 1], my_mesh_2D))

        u_temp = Q[:,2] ./ Q[:,1]
        v_temp = Q[:,3] ./ Q[:,1]
        U_vector = hcat(u_temp, v_temp)
                        
        vector_data = [U_vector] 
        vector_names = ["U"]

        WSE = Q[:,1] + zb_cell
            
        scalar_data = [Q[:,1], Q[:,2], Q[:,3], zb_cell, WSE]
        scalar_names = ["h", "hu", "hv", "zb_cell", "WSE"]

        vtk_fileName = @sprintf("solution_%04d_AdHydraulics.vtk", index)
            
        file_path = joinpath(save_path, vtk_fileName ) 
        export_to_vtk_2D(file_path, my_mesh_2D.nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, scalar_data, scalar_names, vector_data, vector_names)    
    
    end

        
    open(joinpath(save_path, "total_water_volume.csv"), "w") do fo
        println(fo, "total_water_volume")
        for volume in total_water_volume
            println(fo, volume)
        end
    end
end

function export_to_vtk_2D(filename, nodeCoordinates, cellNodesList, cellNodesCount, scalar_data, scalar_names, vector_data, vector_names)
    # Open the file for writing
    open(filename, "w") do file
        # Write the VTK header
        println(file, "# vtk DataFile Version 2.0")
        println(file, "2D Unstructured Mesh")
        println(file, "ASCII")
        println(file, "DATASET UNSTRUCTURED_GRID")
        
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



