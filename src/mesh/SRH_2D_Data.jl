module SRH2D

using LinearAlgebra
using Statistics
using DelimitedFiles

# Class for handling SRH-2D's srhhydro file
struct SRH_2D_SRHHydro
    filename::String
    content::Dict{String, Any}

    function SRH_2D_SRHHydro(filename::String)
        content = Dict{String, Any}()
        obj = new(filename, content)
        parse_srhhydro!(obj)
        return obj
    end
end

function parse_srhhydro!(obj::SRH_2D_SRHHydro)
    content = Dict{String, Any}()
    mannings_n = Dict{Int, Float64}()
    boundary_conditions = Dict{Int, String}()
    iq_params = Dict{Int, Vector{Any}}()

    for line in eachline(obj.filename)
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        parts = split(line, r"\s+")
        if parts[1] == "ManningsN"
            mannings_n[parse(Int, parts[2])] = parse(Float64, parts[3])
        elseif parts[1] == "BC"
            boundary_conditions[parse(Int, parts[2])] = parts[3]
        elseif parts[1] == "IQParams"
            iq_params[parse(Int, parts[2])] = parts[3:end]
        else
            content[parts[1]] = join(parts[2:end], " ")
        end
    end
    content["ManningsN"] = mannings_n
    content["BC"] = boundary_conditions
    content["IQParams"] = iq_params
    obj.content = content
end

# Class for handling SRH-2D's srhgeom file
struct SRH_2D_SRHGeom
    filename::String
    boundary_conditions::Dict{Int, String}
    elements::Vector{Vector{Int}}
    nodes::Matrix{Float64}
    num_elements::Int
    num_nodes::Int

    function SRH_2D_SRHGeom(filename::String, bc::Dict{Int, String})
        elements = Vector{Vector{Int}}()
        nodes = Matrix{Float64}(undef, 0, 3)
        obj = new(filename, bc, elements, nodes, 0, 0)
        parse_geom!(obj)
        return obj
    end
end

function parse_geom!(obj::SRH_2D_SRHGeom)
    nodes = Float64[]
    elements = Vector{Vector{Int}}()
    for line in eachline(obj.filename)
        line = strip(line)
        if startswith(line, "Node")
            parts = split(line)
            append!(nodes, parse.(Float64, parts[2:end]))
        elseif startswith(line, "Elem")
            parts = split(line)
            push!(elements, parse.(Int, parts[2:end]))
        end
    end
    obj.num_nodes = length(nodes) ÷ 3
    obj.num_elements = length(elements)
    obj.nodes = reshape(nodes, 3, obj.num_nodes)'
    obj.elements = elements
end

# Class for handling SRH-2D's srhmat file
struct SRH_2D_SRHMat
    filename::String
    materials::Dict{Int, Vector{Int}}

    function SRH_2D_SRHMat(filename::String)
        materials = Dict{Int, Vector{Int}}()
        obj = new(filename, materials)
        parse_mat!(obj)
        return obj
    end
end

function parse_mat!(obj::SRH_2D_SRHMat)
    current_material = -1
    for line in eachline(obj.filename)
        line = strip(line)
        if startswith(line, "Material")
            parts = split(line)
            current_material = parse(Int, parts[2])
            obj.materials[current_material] = parse.(Int, parts[3:end])
        elseif current_material != -1
            push!(obj.materials[current_material], parse.(Int, split(line))...)
        end
    end
end

# Combined class to integrate the functionalities
struct SRH_2D_Data
    hydro::SRH_2D_SRHHydro
    geom::SRH_2D_SRHGeom
    mat::SRH_2D_SRHMat
end

function SRH_2D_Data(hydro_file::String)
    hydro = SRH_2D_SRHHydro(hydro_file)
    geom_file = hydro.content["Grid"]
    mat_file = hydro.content["HydroMat"]
    geom = SRH_2D_SRHGeom(geom_file, hydro.content["BC"])
    mat = SRH_2D_SRHMat(mat_file)
    return SRH_2D_Data(hydro, geom, mat)
end

end # module
