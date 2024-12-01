function read_srhgeom(file_path)
    name = ""
    grid_unit = ""
    elements = Dict{Int, Vector{Int}}()
    nodes = Dict{Int, Vector{Float64}}()
    boundaries = Dict{String, Vector{Int}}()

    open(file_path, "r") do file
        for line in eachline(file)
            println(line)
            # Parse header information
            if startswith(line, "Name")
                #name = match(r"""Name\s+"(.*)"""", line).captures[1]
                parts = split(line, ' ', limit=2)
                name = parts[2]

            elseif startswith(line, "GridUnit")
                #grid_unit = match(r"""GridUnit\s+"(.*)"""", line).captures[1]
                parts = split(line, ' ', limit=2)
                grid_unit = parts[2]
            
            # Parse elements
            elseif startswith(line, "Elem")
                data = split(line)
                elem_id = parse(Int, data[2])  # Element ID
                node_IDs = parse.(Int, data[3:end])  # Node indices
                elements[elem_id] = node_IDs
            
                # Parse boundary strings (This should preceed the "Nodes" part below)
            elseif startswith(line, "NodeString")
                data = split(line, limit=3)
                println(data)
                #boundary_name = match(r"""NS\s+"(.*)"""", line).captures[1]
                boundary_name = data[2]
                boundary_nodes = parse.(Int, split(data[end]))
                boundaries[boundary_name] = boundary_nodes
            
             # Parse nodes
            elseif startswith(line, "Node")
                data = split(line)
                node_id = parse(Int, data[2])  # Node ID
                coords = parse.(Float64, data[3:5])  # x, y, z
                nodes[node_id] = coords

            end
            
        end
    end

    return name, grid_unit, nodes, elements, boundaries
end

# Function to compute geometric and topological properties
function compute_mesh_properties_srhgeom(numOfCells, numOfFaces, numOfNodes, nodeCoordinates, cellNodesList, cellNodesCount, faceNodes_r_Dict)
    cell_areas = zeros(Float64, numOfCells)
    cell_centroids = zeros(Float64, numOfCells, 2)
    cell_normals = Vector{Vector{Vector{Float64}}}(undef, numOfCells)
    face_normals = Vector{Vector{Float64}}(undef, numOfFaces)
    face_lengths = Vector{Float64}(undef, numOfFaces)   #lenght of the face (edge)

    #loop over all cells
    for cellID in 1:length(cellNodesCount)
        
        #number of nodes for the current cell
        nNodes = cellNodesCount[cellID]

        # Get the node IDs for the current cell
        nodeList = cellNodesList[cellID,:]

        # Get coordinates of the cell's nodes
        vertices = [nodeCoordinates[nodeList[iNode],:][1:2] for iNode in 1:nNodes]

        # Compute area using the shoelace formula
        area = compute_polygon_area_srhgeom(vertices)
        cell_areas[cellID] = area

        # Compute centroid
        centroid = compute_polygon_centroid_srhgeom(vertices)
        cell_centroids[cellID,:] = centroid

        # Compute cell normals (outward)
        normals = compute_face_normals_srhgeom(nNodes, nodeList, vertices)
        cell_normals[cellID] = normals
    end

    #loop over all faces to compute face normals
    for faceID in 1:numOfFaces
        # Get the node IDs for the current face
        nodeList = faceNodes_r_Dict[faceID]

        # Get coordinates of the face's nodes
        vertices = [nodeCoordinates[nodeList[iNode],:][1:2] for iNode in 1:2]

        # Compute edge normals
        face_normal, face_length = compute_face_normal_srhgeom(vertices)
        face_normals[faceID] = face_normal
        face_lengths[faceID] = face_length
    end

    return cell_areas, cell_centroids, cell_normals, face_normals, face_lengths
end

# compute the cell normals
function compute_face_normals_srhgeom(nNodes, nodeList, vertices)
    
    face_normals = Vector{Vector{Float64}}(undef, nNodes)

    for iNode in 1:nNodes
        # Get the current and next node (cyclic)
        jNode = iNode % nNodes + 1
        x1, y1 = vertices[iNode]
        x2, y2 = vertices[jNode]

        # Compute the face vector
        face_vector = [x2 - x1, y2 - y1]
        
        # Compute the outward normal (perpendicular vector)
        normal_vector = [face_vector[2], -face_vector[1]]
        
        # Normalize the normal vector
        normal_vector /= sqrt(normal_vector[1]^2 + normal_vector[2]^2)
        
        # Store the normal vector
        face_normals[iNode] = normal_vector
    end

    return face_normals
    
end

function check_cell_nodes_counter_clockwise_srhgeom(numOfCells, cellNodesList, nodeCoordinates, cellNodesCount)
    #loop over all cells
    for cellID in 1:numOfCells
        #number of nodes for the current cell
        nNodes = cellNodesCount[cellID]
        
        # Get the node IDs for the current cell
        nodeList = cellNodesList[cellID,:]

        # Get coordinates of the cell's nodes
        vertices = [nodeCoordinates[nodeList[iNode],:][1:2] for iNode in 1:nNodes]

        # Compute area using the shoelace formula
        area = compute_polygon_area_srhgeom(vertices)
        if area < 0
            println("Cell $cellID is not counter-clockwise")
            println("Reorering the nodes for cell not implemented yet")
            exit(1)
        end
    end
end

# Helper function: Compute polygon area
function compute_polygon_area_srhgeom(vertices)
    num_vertices = length(vertices)
    area = 0.0
    for i in 1:num_vertices
        x1, y1 = vertices[i]
        x2, y2 = vertices[mod1(i + 1, num_vertices)]
        area += x1 * y2 - x2 * y1
    end

    if area < 0
        println("Area is negative: the polygon is not counter-clockwise")
        println("Reorering the nodes for cell not implemented yet")
        exit(1)
    end

    return abs(area) / 2
end

# Helper function: Compute polygon centroid
function compute_polygon_centroid_srhgeom(vertices)
    num_vertices = length(vertices)
    cx, cy, area_sum = 0.0, 0.0, 0.0
    for i in 1:num_vertices
        x1, y1 = vertices[i]
        x2, y2 = vertices[mod1(i + 1, num_vertices)]
        cross = x1 * y2 - x2 * y1
        cx += (x1 + x2) * cross
        cy += (y1 + y2) * cross
        area_sum += cross
    end
    area = abs(area_sum) / 2
    cx /= 6 * area
    cy /= 6 * area
    return [cx, cy]
end

# Helper function: Compute face normal
function compute_face_normal_srhgeom(vertices)
    num_vertices = length(vertices)

    if num_vertices != 2
        println("Face should have exactly 2 vertices")
        exit(1)
    end

    x1, y1 = vertices[1]
    x2, y2 = vertices[2]
    nx, ny = y2 - y1, -(x2 - x1)  # Perpendicular to the edge
    face_length = sqrt(nx^2 + ny^2)
    face_normal = [nx / face_length, ny / face_length]

    return face_normal, face_length
end

function read_srhmat(file_path)
    num_materials = 0
    material_names = Dict{Int, String}()
    material_elements = Dict{Int, Vector{Int}}()

    open(file_path, "r") do file
        current_material = 0
        for line in eachline(file)
            # Parse number of materials
            if startswith(line, "NMaterials")
                num_materials = parse(Int, split(line)[2])
            
            # Parse material names
            elseif startswith(line, "MatName")
                #data = match(r"""MatName\s+(\d+)\s+"(.*)"""", line)
                data = split(line, limit=3)

                material_id = parse(Int, data[2])
                material_name = data[3]
                material_names[material_id] = material_name

            # Parse material elements
            elseif startswith(line, "Material")
                data = split(line, limit=3)
                current_material = parse(Int, data[2])
                material_elements[current_material] = Int[]

                # Collect element IDs for the first material
                parts = split(data[3])
                element_ids = parse.(Int, parts)
                append!(material_elements[current_material], element_ids)

            # Collect element IDs for the current material
            elseif current_material > 0
                element_ids = parse.(Int, split(line))
                append!(material_elements[current_material], element_ids)
            end
        end
    end

    return num_materials, material_names, material_elements
end


# Function to read and parse SRH-2D hydro file
function read_srhhydro(file_path)
    # Open the file and read lines
    file = open(file_path, "r")
    lines = readlines(file)
    close(file)

    # Initialize a dictionary to store parsed content
    parsed_data = Dict{String, Any}()

    # Iterate through each line and parse based on keywords
    for line in lines
        # Skip empty lines or comments
        if isempty(line) || startswith(line, "#")
            continue
        end

        # Split the line into keyword and value
        split_line = split(line, r"\s+", limit=2)
        if length(split_line) == 2
            keyword, value = split_line[1], split_line[2]
            # Store in the dictionary
            parsed_data[keyword] = value
        elseif length(split_line) == 1
            # For keywords with no explicit value
            parsed_data[split_line[1]] = ""
        end
    end

    return parsed_data
end

