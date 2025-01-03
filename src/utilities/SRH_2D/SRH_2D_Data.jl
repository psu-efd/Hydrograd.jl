"""
A class for SRH-2D data I/O, manipulation, and format conversion
"""
mutable struct SRH_2D_Data
    srhhydro_filename::String
    srhgeom_filename::String
    srhmat_filename::String
    
    # Objects for different file types
    srhhydro_obj::SRH_2D_SRHHydro
    srhgeom_obj::SRH_2D_SRHGeom
    srhmat_obj::SRH_2D_SRHMat
    
    # Manning's n values
    ManningN_cell::Vector{Float64}
    ManningN_node::Vector{Float64}
    
    # XMDF data
    xmdfTimeArray_Nodal::Union{Vector{Float64}, Nothing}
    xmdfAllData_Nodal::Dict{String, Array{Float64}}
    xmdfTimeArray_Cell::Union{Vector{Float64}, Nothing}
    xmdfAllData_Cell::Dict{String, Array{Float64}}

    function SRH_2D_Data(srhhydro_filename::String)
        # Check file existence
        if !isfile(srhhydro_filename)
            throw(ErrorException("The SRHHYDRO file $srhhydro_filename does not exist"))
        end

        # Extract path from srhhydro_filename
        file_path = dirname(srhhydro_filename)

        # Create SRH_2D_SRHHydro object
        srhhydro_obj = SRH_2D_SRHHydro(srhhydro_filename)

        # Get the srhgeom_filename and srhmat_filename and strip any extra quotes
        srhgeom_filename = if isempty(file_path)
            String(strip(srhhydro_get_grid_file_name(srhhydro_obj), ['"', ' ']))
        else
            joinpath(file_path, String(strip(srhhydro_get_grid_file_name(srhhydro_obj), ['"', ' '])))
        end

        srhmat_filename = if isempty(file_path)
            String(strip(srhhydro_get_mat_file_name(srhhydro_obj), ['"', ' ']))
        else
            joinpath(file_path, String(strip(srhhydro_get_mat_file_name(srhhydro_obj), ['"', ' '])))
        end

        @show typeof(srhgeom_filename)
        @show srhgeom_filename

        # Create SRH_2D_SRHGeom and SRH_2D_SRHMat objects
        srhgeom_obj = SRH_2D_SRHGeom(srhgeom_filename, srhhydro_obj.srhhydro_content["BC"])
        srhmat_obj = SRH_2D_SRHMat(srhmat_filename)

        # Initialize Manning's n arrays
        ManningN_cell = zeros(srhgeom_obj.numOfElements)
        ManningN_node = zeros(srhgeom_obj.numOfNodes)

        # Create object
        obj = new(
            srhhydro_filename,
            srhgeom_filename,
            srhmat_filename,
            srhhydro_obj,
            srhgeom_obj,
            srhmat_obj,
            ManningN_cell,
            ManningN_node,
            nothing,  # xmdfTimeArray_Nodal
            Dict{String, Array{Float64}}(),  # xmdfAllData_Nodal
            nothing,  # xmdfTimeArray_Cell
            Dict{String, Array{Float64}}()   # xmdfAllData_Cell
        )

        # Build Manning's n values
        srh_2D_Data_build_manningN_cells_nodes(obj)

        return obj
    end
end

"""
Build Manning's n values at cell centers and nodes
"""
function srh_2D_Data_build_manningN_cells_nodes(data::SRH_2D_Data)
    if gVerbose
        println("Building Manning's n values for cells and nodes in SRH-2D mesh ...")
    end

    # Manning's n dictionary in srhhydro
    nDict = data.srhhydro_obj.srhhydro_content["ManningsN"]
    if gVerbose
        println("nDict = ", nDict)
    end

    # Loop over all cells in the mesh
    for cellI in 1:data.srhgeom_obj.numOfElements
        # Get the material ID of current cell
        matID = srhmat_find_cell_material_ID(data.srhmat_obj, cellI)  # -1 for 0-based indexing
        data.ManningN_cell[cellI] = nDict[matID]
    end

    # Interpolate Manning's n from cell centers to nodes
    for nodeI in 1:data.srhgeom_obj.numOfNodes
        n_temp = 0.0
        
        # Loop over all cells that share the current node
        for i in 1:data.srhgeom_obj.nodeElementsCount[nodeI]
            cellI = data.srhgeom_obj.nodeElementsList[nodeI][i]
            n_temp += data.ManningN_cell[cellI]
        end

        # Take the average
        data.ManningN_node[nodeI] = n_temp / data.srhgeom_obj.nodeElementsCount[nodeI]
    end
end

"""
Get case name from srhhydro file
"""
function srh_2D_Data_get_case_name(data::SRH_2D_Data)
    return srhhydro_get_case_name(data.srhhydro_obj)
end


"""
Read SRH-2D result file in SRHC (cell center) or SRH (point) format
"""
function srh_2D_Data_read_srh_file(data::SRH_2D_Data, srhFileName::String)
    if gVerbose
        println("Reading the SRH/SRHC result file ...")
    end

    # Use DelimitedFiles for CSV reading
    # Note: SRH-2D outputs an extra "," to each line
    df = readdlm(srhFileName, ',', header=true)
    names = df[2][1:end-1]  # Skip last column (empty due to extra comma)
    data_array = df[1]

    return names, data_array
end

"""
Interpolate result from cell center to vertex
"""
function srh_2D_Data_cell_center_to_vertex!(data::SRH_2D_Data, cell_data::Vector{Float64}, vertex_data::Vector{Float64})
    # Initialize vertex data
    fill!(vertex_data, 0.0)
    vertex_counts = zeros(Int, length(vertex_data))

    # Loop through all cells
    for cellI in 1:data.srhgeom_obj.numOfElements
        cell_value = cell_data[cellI]
        
        # Add cell value to all vertices of this cell
        for i in 1:data.srhgeom_obj.elementNodesCount[cellI]
            nodeID = data.srhgeom_obj.elementNodesList[cellI, i]
            vertex_data[nodeID] += cell_value
            vertex_counts[nodeID] += 1
        end
    end

    # Average the values
    for i in eachindex(vertex_data)
        if vertex_counts[i] > 0
            vertex_data[i] /= vertex_counts[i]
        end
    end
end

"""
Interpolate result from vertex to cell center
"""
function srh_2D_Data_vertex_to_cell_center!(data::SRH_2D_Data, vertex_data::Vector{Float64}, cell_data::Vector{Float64})
    # Loop through all cells
    for cellI in 1:data.srhgeom_obj.numOfElements
        cell_value = 0.0
        node_count = data.srhgeom_obj.elementNodesCount[cellI]
        
        # Average values from all vertices of this cell
        for i in 1:node_count
            nodeID = data.srhgeom_obj.elementNodesList[cellI, i]
            cell_value += vertex_data[nodeID]
        end
        
        cell_data[cellI] = cell_value / node_count
    end
end

"""
Get simulation time information
"""
function srh_2D_Data_get_simulation_times(data::SRH_2D_Data)
    return srhhydro_get_simulation_start_end_time(data.srhhydro_obj)
end

"""
Get simulation time step size
"""
function srh_2D_Data_get_simulation_time_step(data::SRH_2D_Data)
    return srhhydro_get_simulation_time_step_size(data.srhhydro_obj)
end

"""
Get number of elements and nodes in the mesh
"""
function srh_2D_Data_get_mesh_size(data::SRH_2D_Data)
    return data.srhgeom_obj.numOfElements, data.srhgeom_obj.numOfNodes
end

"""
Get node coordinates
"""
function srh_2D_Data_get_node_coordinates(data::SRH_2D_Data)
    return data.srhgeom_obj.nodeCoordinates
end

"""
Get element connectivity
"""
function srh_2D_Data_get_element_connectivity(data::SRH_2D_Data)
    return data.srhgeom_obj.elementNodesList, data.srhgeom_obj.elementNodesCount
end


