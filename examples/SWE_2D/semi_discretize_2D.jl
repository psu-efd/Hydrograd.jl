#semin-discretize the 2D SWE to convert it to an ODE system
#This file should be problem specific because we may invert different parameters 

# Problem-specific function to calculate the RHS of ODE: 
# dQdt (= dhdt dq_xdt dq_ydt)
# Q = [h, q_x, q_y] is the solution vector. q_x = h*u_x, q_y = h*u_y
# Q_ghost = [h_ghost, q_x_ghost, q_y_ghost] is the ghost cell values


# In this case, we will invert para = [zb_cells].
function swe_2D_rhs!(dQdt, Q, Q_ghost, para, t, my_mesh_2D, zb_cells, zb_ghostCells, zb_faces, S0,
    swe_2d_constants, ManningN_cells, ManningN_ghostCells,
    nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
    inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length,
    inletQ_TotalA, inletQ_DryWet,  
    nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
    exitH_internalCellIDs, exitH_faceCentroids, exitH_WSE,
    nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, 
    wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals,
    nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs, 
    symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals
    )
    
    #mesh data
    numOfCells = my_mesh_2D.numOfCells
    numOfFaces = my_mesh_2D.numOfFaces
    maxNumOfCellFaces = my_mesh_2D.maxNumOfCellFaces
    
    cell_areas = my_mesh_2D.cell_areas
    cell_normals = my_mesh_2D.cell_normals
    #face_normals = my_mesh_2D.face_normals  
    face_lengths = my_mesh_2D.face_lengths
    
    #faceLeftCellID_Dict = my_mesh_2D.faceLeftCellID_Dict
    #faceRightCellID_Dict = my_mesh_2D.faceRightCellID_Dict

    #get cell neightbors
    cellNeighbors_Dict = my_mesh_2D.cellNeighbors_Dict
    
    g = swe_2d_constants.g
    h_small = swe_2d_constants.h_small
    RiemannSolver = swe_2d_constants.RiemannSolver
    
    h = @view Q[:,1]
    q_x = @view Q[:,2]
    q_y = @view Q[:,3]
    
    #fluxes on faces: 
    fluxes = zeros(Float64, numOfCells, maxNumOfCellFaces, 3)
    flux = zeros(Float64, 3)
    
    #zb at faces and slope at cells 
    #zb_faces = zeros(Float64, numOfFaces)
    #S0 = zeros(Float64, numOfCells, 2)
    
    #println("time t = ", t)
    
    #set the parameter values 
    zb_cells = para

    #interpolate zb from cell to face and compute the bed slope at cells
    interploate_zb_from_cell_to_face_and_compute_S0!(my_mesh_2D, zb_cells, zb_ghostCells, zb_faces, S0)
    
    #Process the boundaries: update ghost cells
    process_inlet_q_boundaries(nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
    inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length,
    inletQ_TotalA, inletQ_DryWet, Q, Q_ghost, swe_2D_constants)
    
    process_exit_h_boundaries(nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
    exitH_internalCellIDs, exitH_faceCentroids, exitH_WSE, Q, Q_ghost, swe_2D_constants)
    
    process_wall_boundaries(nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, 
    wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals, Q, Q_ghost)
    
    process_symmetry_boundaries(nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs, 
    symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals, Q, Q_ghost)
    
    #loop through all celss to calculate the fluxes on faces
    for iCell in 1:numOfCells

        # if iCell==5 && abs(t-5.0)<0.0001
        #     println("iCell = ", iCell)
        #     println("h = ", h[iCell])
        #     println("q_x = ", q_x[iCell])
        #     println("q_y = ", q_y[iCell])
        #     println("S0 = ", S0[iCell,:])
        # end

        #loop through all faces of the cell
        for iFace in 1:my_mesh_2D.cellNodesCount[iCell]
            faceID = abs(my_mesh_2D.cellFacesList[iCell,:][iFace])
            left_cellID = iCell          #left cell is the current cell
            
            right_cellID = abs(cellNeighbors_Dict[iCell][iFace]) #right cell is the neighbor cell (if neighbor cell is negative, it is a ghost cell)

            faceBoundaryID = my_mesh_2D.faceBoundaryID_Dict[faceID]
                        
            #currnet face's normal vector
            face_normal = cell_normals[iCell][iFace]

            if faceBoundaryID==0 #internal face
                hL = h[left_cellID]
                huL = q_x[left_cellID]
                hvL = q_y[left_cellID]
                
                hR = h[right_cellID]
                huR = q_x[right_cellID]
                hvR = q_y[right_cellID]
            else #boundary face
                hL = h[left_cellID]
                huL = q_x[left_cellID]
                hvL = q_y[left_cellID]
                
                hR = Q_ghost[right_cellID, 1]
                huR = Q_ghost[right_cellID, 2]
                hvR = Q_ghost[right_cellID, 3]
            end
            
            if (RiemannSolver == "Roe")
                Riemann_2D_Roe!(flux, hL, huL, hvL, hR, huR, hvR, g, face_normal; hmin=h_small)
                #Riemann_2D_Roe!(flux, h[left_cellID], q_x[left_cellID], q_y[left_cellID], h[right_cellID], q_x[right_cellID], q_y[right_cellID], g, face_normal; hmin=h_small)
                
                #println("Not implemented yet")
                #exit(-1)  #exit with an error code of -1
            elseif (RiemannSolver == "HLL")
                #Riemann_2D_hll!(flux, g, h_face, q_x_face, q_y_face, bcType, h_small)
            elseif (RiemannSolver == "HLLC")
                println("Not implemented yet")
                readline()
                exit(-1)  #exit with an error code of -1
            else
                println("Wrong choice of RiemannSolver")
                readline()
                exit(-1)  #exit with an error code of -1
            end
            
            fluxes[iCell, iFace, 1] = flux[1]
            fluxes[iCell, iFace, 2] = flux[2]
            fluxes[iCell, iFace, 3] = flux[3]

            # if (iCell==5 || iCell==8) && abs(t-5.0)<0.0001
            #     println("iCell = ", iCell, "iFace = ", iFace, "faceID = ", faceID)   
            #     println("h = ", h[iCell])
            #     println("q_x = ", q_x[iCell])
            #     println("q_y = ", q_y[iCell])
            #     println("S0 = ", S0[iCell,:])

            #     println("fluxes = ", fluxes[iCell, iFace, :])
            # end
        end

    end
    
    #initialize the RHS of the ODEs to be zero
    dQdt .*= 0.0
    
    #loop through all cells to calculate the RHS of the ODEs
    for iCell in 1:numOfCells
        
        #get cell area
        cell_area = cell_areas[iCell]
        
        #loop through all faces of the cell
        for iFace in 1:my_mesh_2D.cellNodesCount[iCell]
            faceID = abs(my_mesh_2D.cellFacesList[iCell,:][iFace])
            
            #get face length 
            face_lenght = face_lengths[faceID]
            
            dQdt[iCell,1] -= fluxes[iCell, iFace, 1] * face_lenght 
            dQdt[iCell,2] -= fluxes[iCell, iFace, 2] * face_lenght
            dQdt[iCell,3] -= fluxes[iCell, iFace, 3] * face_lenght
        end
        
        #add flow resistance term and gravity/bed slope term
        if (h[iCell] <= h_small) #if a dry cell, no flow resistance term
            dQdt[iCell,2] += g * h[iCell] * S0[iCell,1] * cell_area
            dQdt[iCell,3] += g * h[iCell] * S0[iCell,2] * cell_area
        else
            u_temp = q_x[iCell] / h[iCell]
            v_temp = q_y[iCell] / h[iCell]
            u_mag = sqrt(u_temp^2 + v_temp^2)

            dQdt[iCell,2] += (g * h[iCell] * S0[iCell,1] - g*ManningN_cells[iCell]^2/h[iCell]^(1.0/3.0)*u_mag*u_temp) * cell_area
            dQdt[iCell,3] += (g * h[iCell] * S0[iCell,2] - g*ManningN_cells[iCell]^2/h[iCell]^(1.0/3.0)*u_mag*v_temp) * cell_area
        end
        
        #divide by cell area
        dQdt[iCell,1] /= cell_area
        dQdt[iCell,2] /= cell_area
        dQdt[iCell,3] /= cell_area
    end 
    
end
