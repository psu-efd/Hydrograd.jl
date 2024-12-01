#2D mesh struct: only non-mutable variables such as node coordinates
mutable struct mesh_2D
    numOfCells::Int64     # Number of cells
    numOfFaces::Int64     # Number of faces 
    numOfNodes::Int64     # Number of nodes
    
    numOfNodeStrings::Int64  # Number of node strings (number of boundaries; it does not include the default boundary)
    
    numOfBoundaries::Int64  # Number of boundaries including the default wall if any.

    numOfTotalCells::Int64  # Total number of cells (including ghost cells)
    
    cellNodesList::Array{Int64,2}  # List of cell's nodes  (2D Int array: [cellID, gMax_Nodes_per_Cell])
    cellNodesCount::Vector{Int64}  # Count of cell's nodes: how many nodes for each cell (1D Int array: [numOfCells])
    
    cellFacesList::Array{Int64,2}  # List of cell's faces  (2D Int array: [cellID, gMax_Nodes_per_Cell]). The numbers of nodes and faces are the same for each cell.
    
    nodeCoordinates::Array{Float64,2}  # Node coordinates: Float64 2D array [numOfNodes, 3]
    
    twoDMeshBoundingbox::Array{Float64,1}  # 2D mesh bounding box: array [xmin, ymin, zmin, xmax, ymax, zmax]
    
    cellBedElevation::Vector{Float64}  # Element bed elevation: Float64 1D array [numOfCells]
    
    nodeStringsDict::Dict{Int64, Vector{Int64}}  #Dictionary of node strings: key: node string ID, value: node string nodes
    
    nodeCellsList::Vector{Vector{Int64}}  # List of node's cells: list of list. Each list contains the cells for each node.
    
    nodeCellsCount::Vector{Int64}  # Count of node's cells: how many cells for each node. 1D Int array: [numOfNodes]
    
    faceNodes_Dict::Dict{Tuple{Int, Int}, Int}  #faceNodes_Dict: key: (node1, node2), value: face ID
    
    faceNodes_r_Dict::Dict{Int, Tuple{Int, Int}}  #faceNodes_r_Dict (reverse of faceNodes_Dict): key: face ID, value: (node1, node2)
    
    faceCells_Dict::Dict{Int, Vector{Int}}  #faceCells_Dict: key: face ID, value: cell ID list 

    faceLeftCellID_Dict::Dict{Int, Int}  #faceLeftCellID_Dict: key: face ID, value: left cell ID

    faceRightCellID_Dict::Dict{Int, Int}  #faceRightCellID_Dict: key: face ID, value: right cell ID

    faceBoundaryID_Dict::Dict{Int, Int}  #faceBoundaryID_Dict: key: face ID, value: boundary ID
    
    boundaryFaces_Dict::Dict{Int, Vector{Int}}  #Dictionary for the List of boundary faces: dictionary: {boundaryID: [list of face IDs]}. It also has the default boundary (default: wall)
    
    allBoundaryFacesIDs_List::Vector{Int}  # List of all boundary faces IDs: all lumped to one list
    
    numOfAllBounaryFaces::Int64  # Number of all boundary faces
    
    boundaryFaces_direction_Dict::Dict{Int, Vector{Int}}  #Dictionary for the direction of boundary faces: {boundaryFaceID: direction}. direction = 1: normal is pointing outward; direction = -1: normal is pointing inward
    
    ghostCellIDs::Vector{Int}  #ghost cell IDs: ghost cells have their own IDs starting from 1
    
    boundaryFaceID_to_ghostCellID_Dict::Dict{Int, Int}  #boundary face ID to ghost cell ID: key: boundary face ID, value: ghost cell ID
    
    boundaryFaceID_to_internalCellID_Dict::Dict{Int, Int}  #boundary face ID to internal cell ID: key: boundary face ID, value: internal cell ID
    
    ghostCellID_to_boundaryFaceID_Dict::Dict{Int, Int}  #ghost cell ID to boundary face ID: key: ghost cell ID, value: boundary face ID
    
    cellNeighbors_Dict::Dict{Int, Vector{Int}}  #cell's neighbors: key: cell ID, value: list of neighbor cell IDs

    cell_areas::Vector{Float64}  #cell areas
    cell_centroids::Array{Float64, 2}  #cell centroids
    cell_normals::Vector{Vector{Vector{Float64}}}  #cell normals
    face_normals::Vector{Vector{Float64}}  #face normals
    face_lengths::Vector{Float64}  #face lengths
    
    
    # Inner constructor with keyword arguments
    function mesh_2D(; 
        numOfCells::Int64=0,
        numOfFaces::Int64=0,
        numOfNodes::Int64=0,
        numOfNodeStrings::Int64=0,
        numOfBoundaries::Int64=0,
        numOfTotalCells::Int64=0,
        cellNodesList::Array{Int64, 2} = Array{Int64}(undef, 0, 0),  # Default to an empty 2D array
        cellNodesCount::Vector{Int64} = Int64[],  # Default to an empty 1D array
        cellFacesList::Array{Int64, 2} = Array{Int64}(undef, 0, 0),  # Default to an empty 2D array
        nodeCoordinates::Array{Float64, 2} = Array{Float64, 2}(undef, 0, 0),  # Default to an uninitialized 3D array with no elements
        twoDMeshBoundingbox::Array{Float64, 1} = Float64[],  # Default to an empty 1D array
        cellBedElevation::Vector{Float64} = Float64[],  # Default to an empty 1D array
        nodeStringsDict::Dict{Int64, Vector{Int64}} = Dict{Int64, Vector{Int64}}(),  # Default to an empty dictionary
        nodeCellsList::Vector{Vector{Int64}} = Vector{Vector{Int64}}(),  # Default to an empty list of lists
        nodeCellsCount::Vector{Int64} = Int64[],  # Default to an empty 1D array
        faceNodes_Dict::Dict{Tuple{Int, Int}, Int} = Dict{Tuple{Int64, Int64}, Int64}(),  # Default to an empty dictionary
        faceNodes_r_Dict::Dict{Int, Tuple{Int, Int}} = Dict{Int, Tuple{Int, Int}}(),  # Default to an empty dictionary
        faceCells_Dict::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}(),  # Default to an empty dictionary
        faceLeftCellID_Dict::Dict{Int, Int} = Dict{Int, Int}(),  # Default to an empty dictionary
        faceRightCellID_Dict::Dict{Int, Int} = Dict{Int, Int}(),  # Default to an empty dictionary
        faceBoundaryID_Dict::Dict{Int, Int} = Dict{Int, Int}(),  # Default to an empty dictionary
        boundaryFaces_Dict::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}(),  # Default to an empty dictionary
        allBoundaryFacesIDs_List::Vector{Int} = Int64[],  # Default to an empty 1D array
        numOfAllBounaryFaces::Int64 = 0,  # Default to 0
        boundaryFaces_direction_Dict::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}(),  # Default to an empty dictionary
        ghostCellIDs::Vector{Int} = Int64[],  # Default to an empty 1D array
        boundaryFaceID_to_ghostCellID_Dict::Dict{Int, Int} = Dict{Int, Int}(),  # Default to an empty dictionary
        boundaryFaceID_to_internalCellID_Dict::Dict{Int, Int} = Dict{Int, Int}(),  # Default to an empty dictionary
        ghostCellID_to_boundaryFaceID_Dict::Dict{Int, Int} = Dict{Int, Int}(),  # Default to an empty dictionary
        cellNeighbors_Dict::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}(),  # Default to an empty dictionary
        cell_areas::Vector{Float64} = Float64[],  # Default to an empty 1D array
        cell_centroids::Array{Float64, 2} = Array{Float64, 2}(undef, 0, 0),  # Default to an uninitialized 2D array with no elements
        cell_normals::Vector{Vector{Vector{Float64}}} = Vector{Vector{Vector{Float64}}}(undef, 0),  # Default to an uninitialized 1D array with no elements
        face_normals::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, 0),  # Default to an uninitialized 1D array with no elements
        face_lengths::Vector{Float64} = Vector{Float64}(undef, 0)  # Default to an uninitialized 1D array with no elements
        )
        
        return new(numOfCells,
                numOfFaces,
                numOfNodes,
                numOfNodeStrings,
                numOfBoundaries,
                numOfTotalCells,
                cellNodesList,
                cellNodesCount,
                cellFacesList,
                nodeCoordinates,
                twoDMeshBoundingbox,
                cellBedElevation,
                nodeStringsDict,
                nodeCellsList,
                nodeCellsCount,
                faceNodes_Dict,
                faceNodes_r_Dict,
                faceCells_Dict,
                faceLeftCellID_Dict,
                faceRightCellID_Dict,
                faceBoundaryID_Dict,
                boundaryFaces_Dict,
                allBoundaryFacesIDs_List,
                numOfAllBounaryFaces,
                boundaryFaces_direction_Dict,
                ghostCellIDs,
                boundaryFaceID_to_ghostCellID_Dict,
                boundaryFaceID_to_internalCellID_Dict,
                ghostCellID_to_boundaryFaceID_Dict,
                cellNeighbors_Dict,
                cell_areas,
                cell_centroids,
                cell_normals,
                face_normals,
                face_lengths
        )
    end
    
end


# Function to initialize the mesh_2D struct
function initialize_mesh_2D(srhgeom_obj, srhhydro_BC) 
    numOfCells = srhgeom_obj.numOfElements    # Number of cells
    numOfNodes = srhgeom_obj.numOfNodes       # Number of nodes
    
    cellNodesList = srhgeom_obj.elementNodesList    # List of cell's nodes  (2D Int array: [cellID, gMax_Nodes_per_Cell])
    cellNodesCount = srhgeom_obj.elementNodesCount  # Count of cell's nodes: how many nodes for each cell (1D Int array: [numOfCells])
    
    cellFacesList = srhgeom_obj.elementEdgesList    # List of cell's faces  (2D Int array: [cellID, gMax_Nodes_per_Cell]). The numbers of nodes and faces are the same for each cell.
    
    println("cellNodesCount: ", cellNodesCount)
    println("cellNodesList: ", cellNodesList)
    println("cellFacesList: ", cellFacesList)
    
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
    
    println("faceNodes_Dict: ")
    for key in keys(faceNodes_Dict)
        println("Key: $key, Value: $(faceNodes_Dict[key])")
    end
    
    println("faceNodes_r_Dict: ")
    for key in sort(collect(keys(faceNodes_r_Dict)))
        println("Key: $key, Value: $(faceNodes_r_Dict[key])")
    end
    
    faceCells_Dict = srhgeom_obj.edgeElements       # Dictionary for the List of face's cells: dictionary: {faceID: [cell list]}
    numOfFaces = length(faceCells_Dict)             # Number of faces
    
    println("numOfFaces: ", numOfFaces)
    
    println("faceCells_Dict: ")
    for key in sort(collect(keys(faceCells_Dict)))
        println("Key: $key, Value: $(faceCells_Dict[key])")
    end
    
    boundaryFaces_Dict = srhgeom_obj.boundaryEdges    # Dictionary for the List of boundary faces: dictionary: {boundaryID: [list of face IDs]}. It also has the default boundary (default: wall)
    allBoundaryFacesIDs_List = srhgeom_obj.allBoundaryEdgeIDs  # List of all boundary faces IDs: all lumped to one list
    numOfAllBounaryFaces = length(allBoundaryFacesIDs_List)  # Number of all boundary faces
    
    println("boundaryFaces_Dict: (before) ", boundaryFaces_Dict)
    
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
    
    println("boundaryFaces_Dict: (after) ", boundaryFaces_Dict)
    println("boundaryFaces_direction_Dict: ", boundaryFaces_direction_Dict)
    println("allBoundaryFacesIDs_List: ", allBoundaryFacesIDs_List)
    
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
    
    println("boundaryFaceID_to_ghostCellID_Dict: ", boundaryFaceID_to_ghostCellID_Dict)
    println("boundaryFaceID_to_internalCellID_Dict: ", boundaryFaceID_to_internalCellID_Dict)
    println("ghostCellID_to_boundaryFaceID_Dict: ", ghostCellID_to_boundaryFaceID_Dict)
    
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
                push!(cellNeighbors_Dict[iCell], -ghostCellID)  #negative ghost cell ID is used to indicate the ghost cell
            end
        end
    end
    
    println("cellNeighbors_Dict: ")
    for key in sort(collect(keys(cellNeighbors_Dict)))
        println("Key: $key, Value: $(cellNeighbors_Dict[key])")
    end
    
    #check cell's nodes counter-clockwise
    check_cell_nodes_counter_clockwise_srhgeom(numOfCells, cellNodesList, nodeCoordinates, cellNodesCount)
    
    #compute mesh properties
    cell_areas, cell_centroids, cell_normals, face_normals, face_lengths = compute_mesh_properties_srhgeom(numOfCells, numOfFaces, numOfNodes, nodeCoordinates, cellNodesList, cellNodesCount, faceNodes_r_Dict)
    
    let counter_temp = 0
        println("cell_areas: ")
        for area in cell_areas
            println(area)
            counter_temp += 1
            if counter_temp == 5
                break
            end
        end
    end 
    
    let counter_temp = 0
        println("cell_centroids: ")
        for i in 1:size(cell_centroids, 1)
            println(cell_centroids[i,:])
            counter_temp += 1
            if counter_temp == 5
                break
            end
        end
    end 
    
    let counter_temp = 0
        println("cell_normals: ")
        for normals in cell_normals
            println(counter_temp+1, ": ", normals)
            counter_temp += 1
            if counter_temp == 5
                #break
            end
        end
    end 
    
    let counter_temp = 0
        println("face_normals: ")
        for normals in face_normals
            println(counter_temp+1, ": ", normals)
            counter_temp += 1
            if counter_temp == 5
                #break
            end
        end
    end 
    
    let counter_temp = 0
        println("face_lengths: ")
        for face_length in face_lengths
            println(counter_temp+1, ": ", face_length)
            counter_temp += 1
            if counter_temp == 5
                #break
            end
        end
    end 

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

        println("faceID: ", iFace, ", faceLeftCellID: ", faceLeftCellID_Dict[iFace], ", faceRightCellID: ", faceRightCellID_Dict[iFace], ", 
            faceBoundaryID: ", faceBoundaryID_Dict[iFace])
    end

    println("faceBoundaryID_Dict: ", faceBoundaryID_Dict)

    #create an empty mesh_2D object
    my_mesh_2D = mesh_2D()

    #modify the mesh_2D object
    my_mesh_2D.numOfCells = numOfCells
    my_mesh_2D.numOfFaces = numOfFaces
    my_mesh_2D.numOfNodes = numOfNodes
    my_mesh_2D.numOfNodeStrings = numOfNodeStrings
    my_mesh_2D.numOfBoundaries = numOfBoundaries
    my_mesh_2D.numOfTotalCells = numOfTotalCells
    my_mesh_2D.cellNodesList = cellNodesList
    my_mesh_2D.cellNodesCount = cellNodesCount
    my_mesh_2D.cellFacesList = cellFacesList
    my_mesh_2D.nodeCoordinates = nodeCoordinates
    my_mesh_2D.twoDMeshBoundingbox = twoDMeshBoundingbox
    my_mesh_2D.cellBedElevation = cellBedElevation
    my_mesh_2D.nodeStringsDict = nodeStringsDict
    my_mesh_2D.nodeCellsList = nodeCellsList
    my_mesh_2D.nodeCellsCount = nodeCellsCount
    my_mesh_2D.faceNodes_Dict = faceNodes_Dict
    my_mesh_2D.faceNodes_r_Dict = faceNodes_r_Dict
    my_mesh_2D.faceCells_Dict = faceCells_Dict
    my_mesh_2D.faceLeftCellID_Dict = faceLeftCellID_Dict
    my_mesh_2D.faceRightCellID_Dict = faceRightCellID_Dict
    my_mesh_2D.faceBoundaryID_Dict = faceBoundaryID_Dict
    my_mesh_2D.boundaryFaces_Dict = boundaryFaces_Dict
    my_mesh_2D.allBoundaryFacesIDs_List = allBoundaryFacesIDs_List
    my_mesh_2D.numOfAllBounaryFaces = numOfAllBounaryFaces
    my_mesh_2D.boundaryFaces_direction_Dict = boundaryFaces_direction_Dict
    my_mesh_2D.ghostCellIDs = ghostCellIDs
    my_mesh_2D.boundaryFaceID_to_ghostCellID_Dict = boundaryFaceID_to_ghostCellID_Dict
    my_mesh_2D.boundaryFaceID_to_internalCellID_Dict = boundaryFaceID_to_internalCellID_Dict
    my_mesh_2D.ghostCellID_to_boundaryFaceID_Dict = ghostCellID_to_boundaryFaceID_Dict
    my_mesh_2D.cellNeighbors_Dict = cellNeighbors_Dict
    my_mesh_2D.cell_areas = cell_areas
    my_mesh_2D.cell_centroids = cell_centroids
    my_mesh_2D.cell_normals = cell_normals
    my_mesh_2D.face_normals = face_normals
    my_mesh_2D.face_lengths = face_lengths
        
    return my_mesh_2D
end

