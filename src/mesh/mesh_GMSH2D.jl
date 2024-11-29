
using LinearAlgebra

# Function to read a Gmsh v2 (Legacy) MSH file
function read_gmsh_v2(filename)
    
    println("Reading Gmsh v2 file: $filename")

    if !isfile(filename)
        error("File not found: $filename")
    end

    nodes = []
    elements = []
    
    open(filename, "r") do file
        while !eof(file)
            line = readline(file)
            
            if startswith(line, "\$Nodes")
                num_nodes = parse(Int, readline(file))
                nodes = zeros(num_nodes, 3) # Store x, y, z coordinates
                for i in 1:num_nodes
                    data = split(readline(file))
                    nodes[i, :] = [parse(Float64, data[2]), parse(Float64, data[3]), parse(Float64, data[4])]
                end
                
            elseif startswith(line, "\$Elements")
                num_elements = parse(Int, readline(file))
                elements = []
                for i in 1:num_elements
                    data = split(readline(file))
                    element_type = parse(Int, data[2])
                    # Collect all element nodes (ignore type for now, generalize later)
                    num_tags = parse(Int, data[3])
                    start_idx = 4 + num_tags
                    element_nodes = parse.(Int, data[start_idx:end])
                    push!(elements, element_nodes)
                end
            end
        end
    end
    
    return nodes, elements
end

# Function to compute polygon area
function compute_polygon_area(vertices)
    num_vertices = size(vertices, 1)
    area = 0.0
    for i in 1:num_vertices
        x1, y1 = vertices[i, :]
        x2, y2 = vertices[mod1(i + 1, num_vertices), :]
        area += x1 * y2 - x2 * y1
    end
    return abs(area) / 2
end

# Function to compute polygon centroid
function compute_polygon_centroid(vertices)
    num_vertices = size(vertices, 1)
    cx, cy, area_sum = 0.0, 0.0, 0.0
    for i in 1:num_vertices
        x1, y1 = vertices[i, :]
        x2, y2 = vertices[mod1(i + 1, num_vertices), :]
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

# Function to compute mesh properties
function compute_mesh_properties(nodes, elements)
    num_elements = length(elements)
    element_areas = zeros(num_elements)
    element_centroids = zeros(num_elements, 2)
    element_normals = []

    for (i, element) in enumerate(elements)
        # Extract vertices of the polygon
        vertices = nodes[element, 1:2]

        # Compute area
        area = compute_polygon_area(vertices)
        element_areas[i] = area

        # Compute centroid
        centroid = compute_polygon_centroid(vertices)
        element_centroids[i, :] = centroid

        # Compute edge normals (for later flux computation)
        local_normals = []
        num_vertices = size(vertices, 1)
        for j in 1:num_vertices
            x1, y1 = vertices[j, :]
            x2, y2 = vertices[mod1(j + 1, num_vertices), :]
            nx, ny = y2 - y1, -(x2 - x1) # Perpendicular to the edge
            length = sqrt(nx^2 + ny^2)
            push!(local_normals, (nx / length, ny / length))
        end
        push!(element_normals, local_normals)
    end

    return element_areas, element_centroids, element_normals
end
