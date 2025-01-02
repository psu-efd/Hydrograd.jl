#process the Manning's n values 

function setup_ManningN(settings, my_mesh_2D, srh_all_Dict)

    # Initialize Manning's n values for cells
    ManningN_cells = zeros(Float64, my_mesh_2D.numOfCells)
    
    srhhydro_ManningsN = srh_all_Dict["srhhydro_ManningsN"]  
    srhmat_numOfMaterials = srh_all_Dict["srhmat_numOfMaterials"]  
    srhmat_matNameList = srh_all_Dict["srhmat_matNameList"]  
    srhmat_matZoneCells = srh_all_Dict["srhmat_matZoneCells"]  
    matID_cells = srh_all_Dict["matID_cells"] 

    if settings.bVerbose
        println("srhhydro_ManningsN: ", srhhydro_ManningsN)
        println("srhmat_numOfMaterials: ", srhmat_numOfMaterials)
        println("srhmat_matNameList: ", srhmat_matNameList)
        println("srhmat_matZoneCells: ", srhmat_matZoneCells)
        println("matID_cells: ", matID_cells)
    end
    
    #loop over all cells to setup Manning's n
    for iCell in 1:my_mesh_2D.numOfCells
        matID = matID_cells[iCell]
        matName = srhmat_matNameList[matID]
        
        if haskey(srhhydro_ManningsN, matID)
            ManningN_cells[iCell] = srhhydro_ManningsN[matID]
        else
            error("Material $matName does not have Manning's n. Assign the default value 0.03.")
        end
    end
    
    #println("ManningN_cells: ", ManningN_cells)

    return ManningN_cells
end

#update Manning's n values based on the provided Manning's n values for each material (zone)
# new_ManningN_values is a vector of Manning's n values for each material (zone)
function update_ManningN(my_mesh_2D::mesh_2D, 
    srh_all_Dict::Dict{String, Any}, 
    new_ManningN_values::Vector{T}) where T <: Real

    matID_cells = srh_all_Dict["matID_cells"]  #material ID for each cell (0-based): 0-default material, 1-first material, 2-second material, etc.

    # Create array directly with comprehension
    ManningN_cells = [new_ManningN_values[matID_cells[i]+1] for i in 1:my_mesh_2D.numOfCells]  #+1 to make matID_cells 1-based (to be consistent with the new_ManningN_values)

    #hack for debugging
    #ManningN_cells = [new_ManningN_values[1] for iCell in 1:my_mesh_2D.numOfCells]
    #ManningN_cells = copy(new_ManningN_values)

    # Check if ManningN_cells and ManningN_ghostCells are AbstractArrays
    #@assert isa(ManningN_cells, AbstractArray) "ManningN_cells must be an AbstractArray"

    # Check if ManningN_cells and ManningN_ghostCells have real number elements
    #@assert eltype(ManningN_cells) <: Real "ManningN_cells must have elements of a real number type"

    return ManningN_cells
end
