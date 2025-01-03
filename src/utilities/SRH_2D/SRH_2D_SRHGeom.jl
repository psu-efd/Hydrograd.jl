# Add these constants if not already defined
const gMax_Nodes_per_Element = 8
const gMax_Elements_per_Node = 20

# VTK cell type mapping
const vtkCellTypeMap = Dict(
    3 => 5,  # triangle
    4 => 9   # quad
)

"""
A class to handle srhgeom file for SRH-2D
"""
mutable struct SRH_2D_SRHGeom
    srhgeom_filename::String
    Name::String
    GridUnit::String
    bcDict::Dict{Int, String}
    
    # Mesh information
    numOfElements::Int
    numOfNodes::Int
    numOfNodeStrings::Int
    
    # Element data
    elementNodesList::Matrix{Int}
    elementNodesCount::Vector{Int}
    elementEdgesList::Matrix{Int}
    vtkCellTypeCode::Vector{Int}
    
    # Node data
    nodeCoordinates::Matrix{Float64}
    twoDMeshBoundingbox::Vector{Float64}
    elementBedElevation::Vector{Float64}
    nodeStringsDict::Dict{Int, Vector{Int}}
    
    # Node-element connectivity
    nodeElementsList::Vector{Vector{Int}}
    nodeElementsCount::Vector{Int}
    
    # Edge data
    edges::Dict{Tuple{Int,Int}, Int}
    edges_r::Dict{Int, Tuple{Int,Int}}
    edgeElements::Dict{Int, Vector{Int}}
    boundaryEdges::Dict{Int, Vector{Int}}
    allBoundaryEdgeIDs::Vector{Int}

    function SRH_2D_SRHGeom(settings, srhgeom_filename::String, bcDict::Dict{Int, String})
        # Initialize with default values
        obj = new(
            srhgeom_filename,
            "", # Name
            "", # GridUnit
            bcDict,
            -1, # numOfElements
            -1, # numOfNodes
            -1, # numOfNodeStrings
            Matrix{Int}(undef, 0, 0), # elementNodesList
            Int[], # elementNodesCount
            Matrix{Int}(undef, 0, 0), # elementEdgesList
            Int[], # vtkCellTypeCode
            Matrix{Float64}(undef, 0, 0), # nodeCoordinates
            Float64[], # twoDMeshBoundingbox
            Float64[], # elementBedElevation
            Dict{Int, Vector{Int}}(), # nodeStringsDict
            Vector{Int}[], # nodeElementsList
            Int[], # nodeElementsCount
            Dict{Tuple{Int,Int}, Int}(), # edges
            Dict{Int, Tuple{Int,Int}}(), # edges_r
            Dict{Int, Vector{Int}}(), # edgeElements
            Dict{Int, Vector{Int}}(), # boundaryEdges
            Int[] # allBoundaryEdgeIDs
        )
        
        # First read to compute number of elements and nodes
        srhgeom_compute_num_of_elements_nodes(settings, obj)
        
        # Initialize arrays with correct sizes
        obj.elementNodesList = zeros(Int, obj.numOfElements, gMax_Nodes_per_Element)
        obj.elementNodesCount = zeros(Int, obj.numOfElements)
        obj.elementEdgesList = zeros(Int, obj.numOfElements, gMax_Nodes_per_Element)
        obj.vtkCellTypeCode = zeros(Int, obj.numOfElements)
        obj.nodeCoordinates = zeros(obj.numOfNodes, 3)
        obj.elementBedElevation = zeros(obj.numOfElements)
        
        # Read the full mesh information
        srhgeom_read_srhgeom_file(settings, obj)
        
        # Build node-element connectivity
        srhgeom_build_node_elements(settings, obj)
        
        # Build edges and boundary edges
        srhgeom_build_edges_and_boundary_edges(settings, obj)
        
        return obj
    end
end

"""
Compute the number of elements and nodes in srhgeom mesh file
"""
function srhgeom_compute_num_of_elements_nodes(settings, geom::SRH_2D_SRHGeom)
    if settings.bVerbose
        println("Computing numbers of elements and nodes from the SRHGEOM file ...")
    end

    elemCount = 0
    nodeCount = 0
    nodeStringCount = 0

    for line in eachline(geom.srhgeom_filename)
        parts = split(strip(line))
        
        if isempty(parts)
            continue
        end

        if parts[1] == "Elem"
            elemCount += 1
        elseif parts[1] == "Node"
            nodeCount += 1
        elseif parts[1] == "NodeString"
            nodeStringCount += 1
        end
    end

    geom.numOfElements = elemCount
    geom.numOfNodes = nodeCount
    geom.numOfNodeStrings = nodeStringCount

    if settings.bVerbose
        println("There are $elemCount elements, $nodeCount nodes, and $nodeStringCount node strings in the mesh.")
    end
end

"""
Get the number of elements and nodes in srhgeom mesh file
"""
function srhgeom_get_num_of_elements_nodes(settings, geom::SRH_2D_SRHGeom)
    if settings.bVerbose
        println("Getting numbers of elements and nodes from the SRHGEOM file ...")
    end

    return geom.numOfElements, geom.numOfNodes
end

"""
Read the SRHGEOM file and populate mesh information
"""
function srhgeom_read_srhgeom_file(settings, geom::SRH_2D_SRHGeom)
    if settings.bVerbose
        println("Reading the SRHGEOM file ...")
    end

    currentNodeStringID = -1

    for line in eachline(geom.srhgeom_filename)
        parts = split(strip(line))
        
        if isempty(parts)
            continue
        end

        if parts[1] == "Elem"
            elemID = parse(Int, parts[2])
            nodes = parse.(Int, parts[3:end])
            geom.elementNodesList[elemID, 1:length(nodes)] = nodes
            geom.elementNodesCount[elemID] = length(nodes)
            geom.vtkCellTypeCode[elemID] = get(vtkCellTypeMap, length(nodes), 0)
        
        elseif parts[1] == "Node"
            nodeID = parse(Int, parts[2])
            coords = parse.(Float64, parts[3:end])
            geom.nodeCoordinates[nodeID, :] = coords
        
        elseif parts[1] == "NodeString"
            nodeStringID = parse(Int, parts[2])
            nodes = parse.(Int, parts[3:end])
            geom.nodeStringsDict[nodeStringID] = nodes
            currentNodeStringID = nodeStringID
        
        elseif lowercase(parts[1]) == "name"
            geom.Name = parts[2]
        
        elseif lowercase(parts[1]) == "gridunit"
            geom.GridUnit = parts[2]
        
        elseif !occursin("srhgeom", lowercase(parts[1]))
            # Assume this is still nodeString
            if currentNodeStringID > 0
                append!(geom.nodeStringsDict[currentNodeStringID], parse.(Int, parts))
            end
        end
    end

    # Calculate bed elevation at cell centers
    for cellI in 1:geom.numOfElements
        elev_temp = 0.0
        for nodeI in 1:geom.elementNodesCount[cellI]
            node_idx = geom.elementNodesList[cellI, nodeI]
            elev_temp += geom.nodeCoordinates[node_idx, 3]
        end
        geom.elementBedElevation[cellI] = elev_temp / geom.elementNodesCount[cellI]
    end

    # Calculate 2D mesh bounding box
    geom.twoDMeshBoundingbox = [
        minimum(geom.nodeCoordinates[:, 1]),  # xmin
        minimum(geom.nodeCoordinates[:, 2]),  # ymin
        minimum(geom.nodeCoordinates[:, 3]),  # zmin
        maximum(geom.nodeCoordinates[:, 1]),  # xmax
        maximum(geom.nodeCoordinates[:, 2]),  # ymax
        maximum(geom.nodeCoordinates[:, 3])   # zmax
    ]

    if settings.bVerbose
        println("2D mesh's bounding box = ", geom.twoDMeshBoundingbox)
    end
end

"""
Build node's element list for all nodes
"""
function srhgeom_build_node_elements(settings, geom::SRH_2D_SRHGeom)
    if settings.bVerbose
        println("Building mesh's node, elements, and topology ...")
    end

    # Initialize nodeElementsList with empty vectors
    geom.nodeElementsList = [Int[] for _ in 1:geom.numOfNodes]
    geom.nodeElementsCount = zeros(Int, geom.numOfNodes)

    # Loop over all cells
    for cellI in 1:geom.numOfElements
        # Loop over all nodes of current cell
        for i in 1:geom.elementNodesCount[cellI]
            nodeID = geom.elementNodesList[cellI, i]
            push!(geom.nodeElementsList[nodeID], cellI)
        end
    end

    # Count elements for each node
    for nodeI in 1:geom.numOfNodes
        geom.nodeElementsCount[nodeI] = length(geom.nodeElementsList[nodeI])
    end
end

"""
Build edges and boundaryEdges dictionaries
"""
function srhgeom_build_edges_and_boundary_edges(settings, geom::SRH_2D_SRHGeom)
    current_edgeID = 1

    # Loop over all elements
    for cellI in 1:geom.numOfElements
        # Loop over all edges of current element
        for i in 1:geom.elementNodesCount[cellI]
            # Get the node IDs for current edge
            nodeID_1 = geom.elementNodesList[cellI, i]
                        
            nodeID_2 = if i != geom.elementNodesCount[cellI]
                geom.elementNodesList[cellI, i+1]
            else
                geom.elementNodesList[cellI, 1]
            end

            if settings.bVerbose
                #println("Cell: $cellI, edge: $i, node IDs: $nodeID_1, $nodeID_2")
            end

            curr_edge_node_IDs = Tuple(sort([nodeID_1, nodeID_2]))

            # Check if edge exists
            if !haskey(geom.edges, curr_edge_node_IDs)
                # Add new edge
                geom.edges[curr_edge_node_IDs] = current_edgeID
                geom.edges_r[current_edgeID] = curr_edge_node_IDs
                geom.edgeElements[current_edgeID] = [cellI]
                current_edgeID += 1
            else
                # Add element to existing edge
                existingEdgeID = geom.edges[curr_edge_node_IDs]
                @assert length(geom.edgeElements[existingEdgeID]) == 1
                push!(geom.edgeElements[existingEdgeID], cellI)
            end
        end
    end

    # Find boundary edges (edges with only one element)
    geom.allBoundaryEdgeIDs = [edge for (edge, elements) in geom.edgeElements if length(elements) == 1]

    if settings.bVerbose
        #println("Boundary edges: ", geom.allBoundaryEdgeIDs)
    end

    # Initialize boundary edge usage flags
    allBoundaryEdgeUsageFlag = Dict(edgeID => false for edgeID in geom.allBoundaryEdgeIDs)

    if settings.bVerbose
        #println("Boundary edge usage flags: ", allBoundaryEdgeUsageFlag)
    end

    # Process boundaries from nodeStrings
    for (nodeString, nodeList) in geom.nodeStringsDict
        #not all NodeStrings are boundaries. Could be monitoring lines. Need to exclude them.
        if !haskey(geom.bcDict, nodeString) || 
           contains(geom.bcDict[nodeString], "WEIR") || 
           contains(geom.bcDict[nodeString], "PRESSURE")
            continue
        end

        current_boundary_edge_list = Int[]

        # Process edges in current boundary
        for i in 1:(length(nodeList)-1)
            nodeID_1 = nodeList[i]
            nodeID_2 = nodeList[i+1]
            
            edge_forward = Tuple([nodeID_1, nodeID_2])
            edge_reverse = Tuple([nodeID_2, nodeID_1])

            if haskey(geom.edges, edge_forward)
                push!(current_boundary_edge_list, geom.edges[edge_forward])
                allBoundaryEdgeUsageFlag[geom.edges[edge_forward]] = true
            elseif haskey(geom.edges, edge_reverse)
                push!(current_boundary_edge_list, -geom.edges[edge_reverse])
                allBoundaryEdgeUsageFlag[geom.edges[edge_reverse]] = true
            else
                error("Boundary edge $edge_forward in NodeString $nodeString not found in edge list. Mesh is wrong.")
            end
        end

        if settings.bVerbose
            #println("Boundary edge list for NodeString $nodeString: ", current_boundary_edge_list)
        end

        geom.boundaryEdges[nodeString] = current_boundary_edge_list
    end

    if settings.bVerbose
        #println("Boundary edges before adding default wall boundary: ", geom.boundaryEdges)
    end

    # Handle unused boundary edges (default wall boundary)
    unusedBoundaryEdges = [edge for edge in geom.allBoundaryEdgeIDs if !allBoundaryEdgeUsageFlag[edge]]
    if !isempty(unusedBoundaryEdges)
        defaultWallBoundaryID = length(geom.boundaryEdges) + 1
        geom.boundaryEdges[defaultWallBoundaryID] = unusedBoundaryEdges
    end

    if settings.bVerbose
        #println("Boundary edges after adding default wall boundary: ", geom.boundaryEdges)
    end

    # Build elementEdgesList
    for cellI in 1:geom.numOfElements
        for i in 1:geom.elementNodesCount[cellI]
            nodeID_1 = geom.elementNodesList[cellI, i]
            
            nodeID_2 = if i != geom.elementNodesCount[cellI]
                geom.elementNodesList[cellI, i+1]
            else
                geom.elementNodesList[cellI, 1]
            end

            edge_forward = Tuple([nodeID_1, nodeID_2])
            edge_reverse = Tuple([nodeID_2, nodeID_1])

            if haskey(geom.edges, edge_forward)
                geom.elementEdgesList[cellI, i] = geom.edges[edge_forward]
            elseif haskey(geom.edges, edge_reverse)
                geom.elementEdgesList[cellI, i] = -geom.edges[edge_reverse]
            else
                error("Cell: $cellI, edge: $i, edge nodes: $edge_forward not found in edge list. Mesh is wrong.")
            end
        end
    end

    if settings.bVerbose
        #println("Element edges list: ", geom.elementEdgesList)
    end
end

