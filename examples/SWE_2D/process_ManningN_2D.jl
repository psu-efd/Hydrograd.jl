#process the Manning's n values 

function setup_ManningN(my_mesh_2D, srh_all_Dict)

    # Initialize Manning's n values for cells
    ManningN_cells = zeros(Float64, my_mesh_2D.numOfCells)
    
    srhhydro_ManningsN = srh_all_Dict["srhhydro_ManningsN"]  
    srhmat_numOfMaterials = srh_all_Dict["srhmat_numOfMaterials"]  
    srhmat_matNameList = srh_all_Dict["srhmat_matNameList"]  
    srhmat_matZoneCells = srh_all_Dict["srhmat_matZoneCells"]  
    matID_cells = srh_all_Dict["matID_cells"] 

    println("srhhydro_ManningsN: ", srhhydro_ManningsN)
    println("srhmat_numOfMaterials: ", srhmat_numOfMaterials)
    println("srhmat_matNameList: ", srhmat_matNameList)
    println("srhmat_matZoneCells: ", srhmat_matZoneCells)
    println("matID_cells: ", matID_cells)
    
    #loop over all cells to setup Manning's n
    for iCell in 1:my_mesh_2D.numOfCells
        matID = matID_cells[iCell]
        matName = srhmat_matNameList[string(matID)]
        
        if haskey(srhhydro_ManningsN, matID)
            ManningN_cells[iCell] = srhhydro_ManningsN[matID]
        else
            error("Material $matName does not have Manning's n. Assign the default value 0.03.")
        end
    end

    #setup Manning's n at ghost cells
    ManningN_ghostCells = update_ghost_cells_scalar(my_mesh_2D, ManningN_cells)    
   
    
    #println("ManningN_cells: ", ManningN_cells)
    #println("ManningN_ghostCells: ", ManningN_ghostCells)

    return ManningN_cells, ManningN_ghostCells
end

#update Manning's n values based on the provided Manning's n values for each material (zone)
# new_ManningN_values is a vector of Manning's n values for each material (zone)
function update_ManningN(my_mesh_2D, srh_all_Dict, new_ManningN_values)

    # Create array directly with comprehension
    ManningN_cells = [new_ManningN_values[matID_cells[i]] for i in 1:my_mesh_2D.numOfCells]

     #update Manning's n at ghost cells
     ManningN_ghostCells = update_ghost_cells_scalar(my_mesh_2D, ManningN_cells)    

     return ManningN_cells, ManningN_ghostCells
end
