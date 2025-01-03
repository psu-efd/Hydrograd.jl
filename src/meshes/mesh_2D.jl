#2D mesh struct: only non-mutable variables 
Base.@kwdef struct mesh_2D
    numOfCells::Int64 = 0     # Number of cells
    numOfFaces::Int64 = 0     # Number of faces 
    numOfNodes::Int64 = 0     # Number of nodes

    maxNumOfCellFaces::Int64 = 0  # Maximum number of faces for each cell
    
    numOfNodeStrings::Int64 = 0  # Number of node strings (number of boundaries; it does not include the default boundary)
    
    numOfBoundaries::Int64 = 0  # Number of boundaries including the default wall if any.

    numOfTotalCells::Int64 = 0  # Total number of cells (including ghost cells)
    
    cellNodesList::Array{Int64,2} = zeros(Int64, 0, 0)  # List of cell's nodes  (2D Int array: [cellID, gMax_Nodes_per_Cell])
    cellNodesCount::Vector{Int64} = zeros(Int64, 0)      # Count of cell's nodes: how many nodes for each cell (1D Int array: [numOfCells])
    
    cellFacesList::Array{Int64,2} = zeros(Int64, 0, 0)  # List of cell's faces  (2D Int array: [cellID, gMax_Nodes_per_Cell]). The numbers of nodes and faces are the same for each cell.
    
    #nodeCoordinates::Array{Float64,2}  # Node coordinates: Float64 2D array [numOfNodes, 3]
    
    twoDMeshBoundingbox::Array{Float64,1} = zeros(Float64, 6)  # 2D mesh bounding box: array [xmin, ymin, zmin, xmax, ymax, zmax]
    
    #cellBedElevation::Vector{Float64}  # Element bed elevation: Float64 1D array [numOfCells]
    
    nodeStringsDict::Dict{Int64, Vector{Int64}} = Dict{Int64, Vector{Int64}}()  #Dictionary of node strings: key: node string ID, value: node string nodes
    
    nodeCellsList::Vector{Vector{Int64}} = Vector{Vector{Int64}}()  # List of node's cells: list of list. Each list contains the cells for each node.
    
    nodeCellsCount::Vector{Int64} = zeros(Int64, 0)  # Count of node's cells: how many cells for each node. 1D Int array: [numOfNodes]
    
    faceNodes_Dict::Dict{Tuple{Int, Int}, Int} = Dict{Tuple{Int, Int}, Int}()  #faceNodes_Dict: key: (node1, node2), value: face ID
    
    faceNodes_r_Dict::Dict{Int, Tuple{Int, Int}} = Dict{Int, Tuple{Int, Int}}()  #faceNodes_r_Dict (reverse of faceNodes_Dict): key: face ID, value: (node1, node2)
    
    faceCells_Dict::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}()  #faceCells_Dict: key: face ID, value: cell ID list 

    faceLeftCellID_Dict::Dict{Int, Int} = Dict{Int, Int}()  #faceLeftCellID_Dict: key: face ID, value: left cell ID

    faceRightCellID_Dict::Dict{Int, Int} = Dict{Int, Int}()  #faceRightCellID_Dict: key: face ID, value: right cell ID

    faceBoundaryID_Dict::Dict{Int, Int} = Dict{Int, Int}()  #faceBoundaryID_Dict: key: face ID, value: boundary ID
    
    boundaryFaces_Dict::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}()  #Dictionary for the List of boundary faces: dictionary: {boundaryID: [list of face IDs]}. It also has the default boundary (default: wall)
    
    allBoundaryFacesIDs_List::Vector{Int} = zeros(Int64, 0)  # List of all boundary faces IDs: all lumped to one list
    
    numOfAllBounaryFaces::Int64 = 0  # Number of all boundary faces
    
    boundaryFaces_direction_Dict::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}()  #Dictionary for the direction of boundary faces: {boundaryFaceID: direction}. direction = 1: normal is pointing outward; direction = -1: normal is pointing inward
    
    ghostCellIDs::Vector{Int} = zeros(Int64, 0)  #ghost cell IDs: ghost cells have their own IDs starting from 1
    
    boundaryFaceID_to_ghostCellID_Dict::Dict{Int, Int} = Dict{Int, Int}()  #boundary face ID to ghost cell ID: key: boundary face ID, value: ghost cell ID
    
    boundaryFaceID_to_internalCellID_Dict::Dict{Int, Int} = Dict{Int, Int}()  #boundary face ID to internal cell ID: key: boundary face ID, value: internal cell ID
    
    ghostCellID_to_boundaryFaceID_Dict::Dict{Int, Int} = Dict{Int, Int}()  #ghost cell ID to boundary face ID: key: ghost cell ID, value: boundary face ID
    
    cellNeighbors_Dict::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}()  #cell's neighbors: key: cell ID, value: list of neighbor cell IDs

    cell_areas::Vector{Float64} = zeros(Float64, 0)  #cell areas
    cell_centroids::Array{Float64, 2} = zeros(Float64, 0, 0)  #cell centroids
    cell_normals::Vector{Vector{Vector{Float64}}} = Vector{Vector{Vector{Float64}}}()  #cell normals
    face_normals::Vector{Vector{Float64}} = Vector{Vector{Float64}}()  #face normals
    face_lengths::Vector{Float64} = zeros(Float64, 0)  #face lengths
    
end


# Function to initialize the mesh_2D struct
function initialize_mesh_2D(srhgeom_obj, srhhydro_BC) 
    numOfCells = srhgeom_obj.numOfElements    # Number of cells
    numOfNodes = srhgeom_obj.numOfNodes       # Number of nodes
    
    cellNodesList = srhgeom_obj.elementNodesList    # List of cell's nodes  (2D Int array: [cellID, gMax_Nodes_per_Cell])
    cellNodesCount = srhgeom_obj.elementNodesCount  # Count of cell's nodes: how many nodes for each cell (1D Int array: [numOfCells])

    maxNumOfCellFaces = maximum(cellNodesCount)  # Maximum number of faces for each cell
    
    cellFacesList = srhgeom_obj.elementEdgesList    # List of cell's faces  (2D Int array: [cellID, gMax_Nodes_per_Cell]). The numbers of nodes and faces are the same for each cell.
    
    # println("cellNodesCount: ", cellNodesCount)
    # println("cellNodesList: ", cellNodesList)
    # println("cellFacesList: ", cellFacesList)
    # println("maxNumOfCellFaces: ", maxNumOfCellFaces)
    
    nodeCoordinates = srhgeom_obj.nodeCoordinates      # Node coordinates: Float64 2D array [numOfNodes, 3]
    twoDMeshBoundingbox = srhgeom_obj.twoDMeshBoundingbox  # 2D mesh bounding box: array [xmin, ymin, zmin, xmax, ymax, zmax]
    cellBedElevation = srhgeom_obj.elementBedElevation  # Element bed elevation: Float64 1D array [numOfCells]
    nodeStringsDict = srhgeom_obj.nodeStringsDict        #Dictionary of node strings: key: node string ID, value: node string nodes
    nodeCellsList = srhgeom_obj.nodeElementsList      # List of node's cells: list of list. Each list contains the cells for each node. 
    nodeCellsCount = srhgeom_obj.nodeElementsCount    # Count of node's cells: how many cells for each node. 1D Int array: [numOfNodes]
    
    # Convert srhgeom_obj.edges to Julia dictionary
    #faceNodes_Dict: key: (node1, node2), value: face ID
    faceNodes_Dict = Dict{Tuple{Int, Int}, Int}(tuple(convert(Int64,k[1]), convert(Int64,k[2])) => convert(Int64,v) for (k, v) in srhgeom_obj.edges)
    #faceNodes_r_Dict (reverse of faceNodes_Dict): key: face ID, value: (node1, node2)
    faceNodes_r_Dict = Dict{Int, Tuple{Int, Int}}(convert(Int64,k) => tuple(convert(Int64,v[1]), convert(Int64,v[2])) for (k, v) in srhgeom_obj.edges_r)
    
    # println("faceNodes_Dict: ")
    # for key in keys(faceNodes_Dict)
    #     println("Key: $key, Value: $(faceNodes_Dict[key])")
    # end
    
    # println("faceNodes_r_Dict: ")
    # for key in sort(collect(keys(faceNodes_r_Dict)))
    #     println("Key: $key, Value: $(faceNodes_r_Dict[key])")
    # end
    
    faceCells_Dict = srhgeom_obj.edgeElements       # Dictionary for the List of face's cells: dictionary: {faceID: [cell list]}
    numOfFaces = length(faceCells_Dict)             # Number of faces
    
    #println("numOfFaces: ", numOfFaces)
    
    #println("faceCells_Dict: ")
    #for key in sort(collect(keys(faceCells_Dict)))
    #    println("Key: $key, Value: $(faceCells_Dict[key])")
    #end
    
    boundaryFaces_Dict = srhgeom_obj.boundaryEdges    # Dictionary for the List of boundary faces: dictionary: {boundaryID: [list of face IDs]}. It also has the default boundary (default: wall)
    allBoundaryFacesIDs_List = srhgeom_obj.allBoundaryEdgeIDs  # List of all boundary faces IDs: all lumped to one list
    numOfAllBounaryFaces = length(allBoundaryFacesIDs_List)  # Number of all boundary faces
    
    #println("boundaryFaces_Dict: (before) ", boundaryFaces_Dict)
    
    #boundaryFaces_Dict's face ID is negative for the boundary faces whose normal is pointing inwward the domain. Here is to record that.
    boundaryFaces_direction_Dict = Dict()  #Dictionary for the direction of boundary faces: {boundaryFaceID: direction}. direction = 1: normal is pointing outward; direction = -1: normal is pointing inward
    
    for key in keys(boundaryFaces_Dict)
        faceIDs = boundaryFaces_Dict[key]
        
        abs_faceIDs = abs.(faceIDs)
        boundaryFaces_Dict[key] = abs_faceIDs                #now boundaryFaces_Dict[key] has all positive face IDs
    end

    #also loop through allBoundaryFacesIDs_List to check if the face ID is negative (a little redundant)
    for (boundaryID, boundaryFace_IDs) in boundaryFaces_Dict
        boundaryFaces_direction = [1 for i in 1:length(boundaryFace_IDs)]  #initialize the direction of boundary faces to 1

        for (iBoundaryFace, boundaryFace_ID) in  enumerate(boundaryFace_IDs)

            if (boundaryFace_ID in cellFacesList) && (-boundaryFace_ID) in cellFacesList    #if the positive and negative face ID are in cellFacesList, it is wrong. A boundary face should only show up once. 
                println("boundary face ID is both positive and negative. It is wrong.")
                readline()
                exit(-1)
            end

            if (-boundaryFace_ID) in cellFacesList    #if the negative face ID is in cellFacesList, it is an inward normal face
                boundaryFaces_direction[iBoundaryFace] = -1
            end
        end

        boundaryFaces_direction_Dict[boundaryID] = boundaryFaces_direction
    end

    #above is the only place the sign of cellFacesList is used. Here, the sign is turned to positive for all. 
    cellFacesList = abs.(cellFacesList)

    #println("cellFacesList after abs: ", cellFacesList)

    
    #println("boundaryFaces_Dict: (after) ", boundaryFaces_Dict)
    #println("boundaryFaces_direction_Dict: ", boundaryFaces_direction_Dict)
    #println("allBoundaryFacesIDs_List: ", allBoundaryFacesIDs_List)
    
    numOfNodeStrings = srhgeom_obj.numOfNodeStrings    # Number of node strings (number of boundaries; it does not include the default boundary)
    
    #srhhydro_BC may not include the defautl wall boundary. Add the default wall boundary to srhhydro_BC if it does not exist.
    if length(srhhydro_BC) < length(boundaryFaces_Dict)
        @assert length(srhhydro_BC) == length(boundaryFaces_Dict) - 1
        
        for key in keys(boundaryFaces_Dict)
            if !haskey(srhhydro_BC, key)
                srhhydro_BC[key] = "wall"
            end
        end
    end
    
    numOfBoundaries = length(srhhydro_BC)    # Number of boundaries including the default wall if any.
    
    #ghost cells for boundary faces: ghost cell's index is the same as their index in allBoundaryFacesIDs_List, e.g., the first ghost cell is for the first boundary face in allBoundaryFacesIDs_List.
    #boundaryFacesGhostCells_Dict = Dict{Int, Tuple{Int, Int}}()  #Dictionary for the ghost cells of boundary faces: {boundaryFaceID: (left ghost cell, right ghost cell)}
    
    #total number of cells (including ghost cells)
    numOfTotalCells = numOfCells + numOfAllBounaryFaces
    
    #ghost cell IDs: ghost cells have their own IDs starting from 1
    ghostCellIDs = [i for i in 1:numOfAllBounaryFaces]
    
    #boundary face ID to ghost cell ID: key: boundary face ID, value: ghost cell ID
    boundaryFaceID_to_ghostCellID_Dict = Dict{Int, Int}()
    #loop over all boundary faces to find the ghost cell ID for each boundary face
    for iBoundaryFace in 1:numOfAllBounaryFaces
        boundaryFaceID_to_ghostCellID_Dict[allBoundaryFacesIDs_List[iBoundaryFace]] = ghostCellIDs[iBoundaryFace]
    end
    
    #boundary face ID to internal cell ID: key: boundary face ID, value: internal cell ID
    boundaryFaceID_to_internalCellID_Dict = Dict{Int, Int}()
    #loop over all boundary faces to find the internal cell ID for each boundary face
    for iBoundaryFace in 1:numOfAllBounaryFaces
        faceID = allBoundaryFacesIDs_List[iBoundaryFace]
        
        faceCells = faceCells_Dict[abs(faceID)]
        
        @assert length(faceCells) == 1  #boundary face has only one cell
        
        boundaryFaceID_to_internalCellID_Dict[allBoundaryFacesIDs_List[iBoundaryFace]] = faceCells[1]
    end
    
    #ghost cell ID to boundary face ID: key: ghost cell ID, value: boundary face ID
    ghostCellID_to_boundaryFaceID_Dict = Dict{Int, Int}()
    #loop over all ghost cells to find the boundary face ID for each ghost cell
    for iGhostCell in 1:numOfAllBounaryFaces
        ghostCellID_to_boundaryFaceID_Dict[ghostCellIDs[iGhostCell]] = allBoundaryFacesIDs_List[iGhostCell]
    end
    
    #println("boundaryFaceID_to_ghostCellID_Dict: ", boundaryFaceID_to_ghostCellID_Dict)
    #println("boundaryFaceID_to_internalCellID_Dict: ", boundaryFaceID_to_internalCellID_Dict)
    #println("ghostCellID_to_boundaryFaceID_Dict: ", ghostCellID_to_boundaryFaceID_Dict)
    
    #cell's neighbors: key: cell ID, value: list of neighbor cell IDs
    cellNeighbors_Dict = Dict{Int, Vector{Int}}()
    
    #loop over all cells to find the neighbors of each cell
    for iCell in 1:numOfCells
        cellNeighbors_Dict[iCell] = []
        
        nNodes = cellNodesCount[iCell]
        
        #loop over all the faces of the current cell
        for iFace in 1:nNodes
            faceID = cellFacesList[iCell, iFace]
            
            faceCells = faceCells_Dict[abs(faceID)]
            
            if length(faceCells) == 2   #internal face
                if faceCells[1] == iCell
                    push!(cellNeighbors_Dict[iCell], faceCells[2])
                else
                    push!(cellNeighbors_Dict[iCell], faceCells[1])
                end
            else  #boundary face
                ghostCellID = boundaryFaceID_to_ghostCellID_Dict[abs(faceID)]
                push!(cellNeighbors_Dict[iCell], ghostCellID)  
            end
        end
    end
    
    #println("cellNeighbors_Dict: ")
    #for key in sort(collect(keys(cellNeighbors_Dict)))
    #    println("Key: $key, Value: $(cellNeighbors_Dict[key])")
    #end
    
    #check cell's nodes counter-clockwise
    check_cell_nodes_counter_clockwise_srhgeom(numOfCells, cellNodesList, nodeCoordinates, cellNodesCount)
    
    #compute mesh properties
    cell_areas, cell_centroids, cell_normals, face_normals, face_lengths = compute_mesh_properties_srhgeom(numOfCells, numOfFaces, numOfNodes, nodeCoordinates, cellNodesList, cellNodesCount, faceNodes_r_Dict)
    
    #let counter_temp = 0
    #    println("cell_areas: ")
    #    for area in cell_areas
    #        println(area)
    #        counter_temp += 1
    #        if counter_temp == 5
    #            break
    #        end
    #    end
    #end 
    
    # let counter_temp = 0
    #     println("cell_centroids: ")
    #     for i in 1:size(cell_centroids, 1)
    #         println(cell_centroids[i,:])
    #         counter_temp += 1
    #         if counter_temp == 5
    #             break
    #         end
    #     end
    # end 
    
    # let counter_temp = 0
    #     println("cell_normals: ")
    #     for normals in cell_normals
    #         println(counter_temp+1, ": ", normals)
    #         counter_temp += 1
    #         if counter_temp == 5
    #             #break
    #         end
    #     end
    # end 
    
    # let counter_temp = 0
    #     println("face_normals: ")
    #     for normals in face_normals
    #         println(counter_temp+1, ": ", normals)
    #         counter_temp += 1
    #         if counter_temp == 5
    #             #break
    #         end
    #     end
    # end 
    
    # let counter_temp = 0
    #     println("face_lengths: ")
    #     for face_length in face_lengths
    #         println(counter_temp+1, ": ", face_length)
    #         counter_temp += 1
    #         if counter_temp == 5
    #             #break
    #         end
    #     end
    # end 

    #compute the left and right cells of each face
    #For internal face, the left cell -> the right cell is in the same direction of the outer normal vector 
    #For boundary face, the left cell is the internal cell and the right cell is the ghost cell becuase 
    #the outer normal vector of a boundary face is pointing outward the domain.
    #loop over all faces
    faceLeftCellID_Dict = Dict{Int, Int}()
    faceRightCellID_Dict = Dict{Int, Int}()
    faceBoundaryID_Dict = Dict{Int, Int}()

    for iFace in keys(faceCells_Dict)
        
        faceCells = faceCells_Dict[iFace]
        
        if length(faceCells) == 2   #internal face
            #get the face normals
            face_normal = face_normals[iFace]

            #get the cell centroids
            cell_centroid1 = cell_centroids[faceCells[1], :]
            cell_centroid2 = cell_centroids[faceCells[2], :]

            #get the vector from cell1 to cell2
            vector_cell1_to_cell2 = (cell_centroid2 - cell_centroid1)[1:2]

            #compute the inner product of the face normal and the vector from cell1 to cell2
            inner_product = dot(face_normal, vector_cell1_to_cell2)

            if inner_product > 0.0
                faceLeftCellID_Dict[iFace] = faceCells[1]
                faceRightCellID_Dict[iFace] = faceCells[2]
            else
                faceLeftCellID_Dict[iFace] = faceCells[2]
                faceRightCellID_Dict[iFace] = faceCells[1]
            end

            faceBoundaryID_Dict[iFace] = 0  #0 means internal face

        else  #boundary face
            #get the ghost cell ID
            ghostCellID = boundaryFaceID_to_ghostCellID_Dict[iFace]

            faceLeftCellID_Dict[iFace] = faceCells[1]
            faceRightCellID_Dict[iFace] = ghostCellID

            #find which boundary the boundary face belongs to
            for (boundaryID, boundaryFaceIDs) in boundaryFaces_Dict
                if iFace in boundaryFaceIDs
                    faceBoundaryID_Dict[iFace] = boundaryID
                    break
                end
            end

            #check whether the boundary face has found its boundary
            if !haskey(faceBoundaryID_Dict, iFace)
                println("Error: boundary face has not found its boundary: ", iFace)
                readline()
                exit(-1)
            end
        end

        #println("faceID: ", iFace, ", faceLeftCellID: ", faceLeftCellID_Dict[iFace], ", faceRightCellID: ", faceRightCellID_Dict[iFace], ", 
        #    faceBoundaryID: ", faceBoundaryID_Dict[iFace])
    end

    #println("faceBoundaryID_Dict: ", faceBoundaryID_Dict)

    #create an empty mesh_2D object
    my_mesh_2D = mesh_2D(
        numOfCells = numOfCells,
        numOfFaces = numOfFaces,
        numOfNodes = numOfNodes,
        maxNumOfCellFaces = maxNumOfCellFaces,
        numOfNodeStrings = numOfNodeStrings,
        numOfBoundaries = numOfBoundaries,
        numOfTotalCells = numOfTotalCells,
        cellNodesList = cellNodesList,
        cellNodesCount = cellNodesCount,
        cellFacesList = cellFacesList,
        twoDMeshBoundingbox = twoDMeshBoundingbox,
        nodeStringsDict = nodeStringsDict,
        nodeCellsList = nodeCellsList,
        nodeCellsCount = nodeCellsCount,
        faceNodes_Dict = faceNodes_Dict,
        faceNodes_r_Dict = faceNodes_r_Dict,
        faceCells_Dict = faceCells_Dict,
        faceLeftCellID_Dict = faceLeftCellID_Dict,
        faceRightCellID_Dict = faceRightCellID_Dict,
        faceBoundaryID_Dict = faceBoundaryID_Dict,
        boundaryFaces_Dict = boundaryFaces_Dict,
        allBoundaryFacesIDs_List = allBoundaryFacesIDs_List,
        numOfAllBounaryFaces = numOfAllBounaryFaces,
        boundaryFaces_direction_Dict = boundaryFaces_direction_Dict,
        ghostCellIDs = ghostCellIDs,
        boundaryFaceID_to_ghostCellID_Dict = boundaryFaceID_to_ghostCellID_Dict,
        boundaryFaceID_to_internalCellID_Dict = boundaryFaceID_to_internalCellID_Dict,
        ghostCellID_to_boundaryFaceID_Dict = ghostCellID_to_boundaryFaceID_Dict,
        cellNeighbors_Dict = cellNeighbors_Dict,
        cell_areas = cell_areas,
        cell_centroids = cell_centroids,
        cell_normals = cell_normals,
        face_normals = face_normals,
        face_lengths = face_lengths
    )
        
    return my_mesh_2D
end

# Function to compute geometric and topological properties
function compute_mesh_properties_srhgeom(numOfCells, numOfFaces, numOfNodes, nodeCoordinates, cellNodesList, cellNodesCount, faceNodes_r_Dict)
    cell_areas = zeros(Float64, numOfCells)
    cell_centroids = zeros(Float64, numOfCells, 2)
    cell_normals = Vector{Vector{Vector{Float64}}}(undef, numOfCells)
    face_normals = Vector{Vector{Float64}}(undef, numOfFaces)
    face_lengths = Vector{Float64}(undef, numOfFaces)   #lenght of the face (edge)

    #loop over all cells
    for cellID in eachindex(cellNodesCount)
        
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

