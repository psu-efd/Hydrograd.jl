"""
A class to handle srhmat file for SRH-2D
"""
mutable struct SRH_2D_SRHMat
    srhmat_filename::String
    numOfMaterials::Int
    matNameList::Dict{String, String}  # ID and name
    matZoneCells::Dict{Int, Vector{Int}}  # material zone ID and cell list

    function SRH_2D_SRHMat(srhmat_filename::String)
        obj = new(
            srhmat_filename,
            -1,
            Dict{String, String}(),
            Dict{Int, Vector{Int}}()
        )
        
        # Add default material name
        obj.matNameList["-1"] = "Default"
        
        # Build material zones data
        srhmat_build_material_zones_data(obj)
        return obj
    end
end

"""
Build the data for material zones from reading the srhmat file
"""
function srhmat_build_material_zones_data(mat::SRH_2D_SRHMat)
    if gVerbose
        println("Reading the SRHMAT file ...")
    end

    # Debug prints
    @show mat.srhmat_filename
    @show isfile(mat.srhmat_filename)
    @show pwd()  # Show current working directory

    # Check file existence before opening
    if !isfile(mat.srhmat_filename)
        throw(ErrorException("SRHMAT file $(mat.srhmat_filename) does not exist. Current directory: $(pwd())"))
    end

    # Initialize result dictionaries
    res_MatNames = Dict{String, String}("-1" => "Default")
    res_Materials = Dict{Int, Vector{Int}}()

    current_MaterialID = -1
    current_MaterialCellList = Int[]

    println("About to open file...")  # Debug print
    
    open(mat.srhmat_filename, "r") do file
        println("File opened successfully")  # Debug print
        println("First line: ", readline(file))  # Try to read first line
        
        # Reset file position to start
        seekstart(file)
        
        for line in eachline(file)
            println("Reading line: ", line)  # Debug print
            parts = split(strip(line))
            
            if isempty(parts)
                continue
            end

            if parts[1] == "SRHMAT"
                continue
            elseif parts[1] == "NMaterials"
                mat.numOfMaterials = parse(Int, parts[2])
            elseif parts[1] == "MatName"
                res_MatNames[parts[2]] = parts[3]
            elseif parts[1] == "Material"
                # Save previous material zone if exists
                if current_MaterialID != -1 && !isempty(current_MaterialCellList)
                    res_Materials[current_MaterialID] = copy(current_MaterialCellList)
                    current_MaterialCellList = Int[]
                end
                
                # Start new material zone
                current_MaterialID = parse(Int, parts[2])
                append!(current_MaterialCellList, parse.(Int, parts[3:end]))
            else
                # Assume continuation of material cell list
                append!(current_MaterialCellList, parse.(Int, parts))
            end
        end
    end

    # Add the last material zone
    if current_MaterialID != -1
        res_Materials[current_MaterialID] = copy(current_MaterialCellList)
    end

    # Add default material with empty cell list
    res_Materials[-1] = Int[]

    if gVerbose
        println("res_MatNames: ", res_MatNames)
        println("res_Materials: ", res_Materials)
    end

    mat.matNameList = res_MatNames
    mat.matZoneCells = res_Materials
end

"""
Find material (Manning's n) ID for a given cell ID
"""
function srhmat_find_cell_material_ID(mat::SRH_2D_SRHMat, cellID::Int)
    # Check if material zones data exists
    if isempty(mat.matZoneCells)
        srhmat_build_material_zones_data(mat)
    end

    # Search for cellID in material zones
    for (matID, cellList) in mat.matZoneCells
        if (cellID + 1) in cellList  # cellID+1 
            return matID
        end
    end

    # If cell not found, warn and return default
    println("pyHMT2D: In find_cell_material_ID(cellID), cellID = $cellID is not found. ",
            "Check mesh and material coverage. Default is used.")
    return 0
end


