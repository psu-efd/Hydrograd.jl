#process the Manning's n values 

function setup_ManningN!(ManningN_cells, ManningN_ghostCells, my_mesh_2D, srh_all_Dict)
    
    srhhydro_ManningsN = srh_all_Dict["srhhydro_ManningsN"]  
    srhmat_numOfMaterials = srh_all_Dict["srhmat_numOfMaterials"]  
    srhmat_matNameList = srh_all_Dict["srhmat_matNameList"]  
    srhmat_matZoneCells = srh_all_Dict["srhmat_matZoneCells"]  
    matID_cells = srh_all_Dict["matID_cells"] 
    
    #loop over all cells to setup Manning's n
    for iCell in 1:my_mesh_2D.numOfCells
        matID = matID_cells[iCell]
        matName = srhmat_matNameList[string(matID)]
        
        if haskey(srhhydro_ManningsN, matID)
            ManningN_cells[iCell] = srhhydro_ManningsN[matID]
        else
            println("Material ", matName, " does not have Manning's n. Assign the default value 0.03.")
            exit(-1)
        end
    end

    #setup Manning's n at ghost cells
    update_ghost_cells_scalar!(my_mesh_2D.numOfAllBounaryFaces, my_mesh_2D.allBoundaryFacesIDs_List, 
        my_mesh_2D.faceCells_Dict, ManningN_cells, ManningN_ghostCells)
   
    #println("ManningN_cells: ", ManningN_cells)
    #println("ManningN_ghostCells: ", ManningN_ghostCells)
end