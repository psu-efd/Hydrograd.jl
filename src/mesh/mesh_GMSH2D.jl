


struct MeshNode
    id::Int
    x::Float64
    y::Float64
    z::Float64
end

struct MeshElement
    id::Int
    type::Int
    tagCount::Int
    nodeIDs::Vector{Int}
end

struct MeshData
    nodes::Vector{MeshNode}
    elements::Vector{MeshElement}
end

"""
    readGmshMesh(filename::String) -> MeshData

Read a Gmsh 2D mesh file and return the mesh data.
"""
function readGmshMesh(filename::String)
    if !isfile(filename)
        error("File not found: $filename")
    end

    nodes = MeshNode[]
    elements = MeshElement[]

    open(filename, "r") do file
        readNodes = false
        readElements = false

        for line in eachline(file)
            if starts_with(line, "$Nodes")
                readNodes = true
                continue
            elseif starts_with(line, "$EndNodes")
                readNodes = false
                continue
            elseif starts_with(line, "$Elements")
                readElements = true
                continue
            elseif starts_with(line, "$EndElements")
                readElements = false
                continue
            end

            if readNodes
                nodeData = split(line)
                push!(nodes, MeshNode(parse(Int, nodeData[1]), parse(Float64, nodeData[2]), parse(Float64, nodeData[3]), parse(Float64, nodeData[4])))
            elseif readElements
                elementData = split(line)
                elementID = parse(Int, elementData[1])
                elementType = parse(Int, elementData[2])
                tagCount = parse(Int, elementData[3])
                nodeIDs = [parse(Int, id) for id in elementData[4:end]]
                push!(elements, MeshElement(elementID, elementType, tagCount, nodeIDs))
            end
        end
    end

    return MeshData(nodes, elements)
end

