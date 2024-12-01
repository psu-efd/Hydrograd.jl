# process initial condition, such as setting up initial free surface, water depth, etc. 
# This file should be problem specific because each problem should have different ICs. 

# setup initial condition: free surface 
function setup_initial_eta!(numOfCells, nodeCoordinates, cellNodesList, cellNodesCount, cell_centroids, eta, zb_cell, h, bPlot::Bool=false)
    # parameters for initial free surface profile setup
    
    xMid = (minimum(nodeCoordinates[:,1])+maximum(nodeCoordinates[:,1])) / 2.0   #mid point of the domain
    
    bump_center_x = xMid  # center of the bump
    
    h_small = 0.01
    
    #loop over cells
    @inbounds for i in 1:numOfCells
        if cell_centroids[i,1] < bump_center_x
            eta[i] = 1.0
        else
            eta[i] = 1.0    #0.5
        end
    end
    
    #update water depth
    @inbounds for i in 1:numOfCells
        h[i] = eta[i] - zb_cell[i]
    end
    
    h[h.<0.0] .= h_small  #ensure positivity of water depth h
    
    #update the free surface elevation again in case h has been clipped
    @inbounds for i in 1:numOfCells
        eta[i] = h[i] + zb_cell[i]
    end
    
    #optionally plot the free surface for checking 
    if bPlot
        vector_data = [] 
        vector_names = []
        
        scalar_data = [eta, h, zb_cell]
        scalar_names = ["eta", "h", "zb_cell"]
        
        file_path = joinpath(@__DIR__, "eta_h_zb.vtk" ) 
        export_to_vtk_2D(file_path, nodeCoordinates, cellNodesList, cellNodesCount, scalar_data, scalar_names, vector_data, vector_names)    
        println("eta, h, and zb are saved to ", file_path)
        #exit(0)
    end
    
end

#update ghost cells for eta, h, q_x, q_y
function update_ghost_cells_eta_h_q!(numOfAllBounaryFaces, allBoundaryFacesIDs_List, faceCells_Dict, eta, h, q_x, q_y, 
    eta_ghostCells, h_ghostCells, q_x_ghostCells, q_y_ghostCells)

    for iBoundaryFace in 1:numOfAllBounaryFaces
        cellID_neighbor = faceCells_Dict[allBoundaryFacesIDs_List[iBoundaryFace]][1]
        
        eta_ghostCells[iBoundaryFace] = eta[cellID_neighbor]
        h_ghostCells[iBoundaryFace] = h[cellID_neighbor]
        q_x_ghostCells[iBoundaryFace] = q_x[cellID_neighbor]
        q_y_ghostCells[iBoundaryFace] = q_y[cellID_neighbor]
    end
    
    println("h_ghostCells: ", h_ghostCells)
end 


