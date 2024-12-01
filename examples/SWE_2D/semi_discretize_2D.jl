#semin-discretize the 2D SWE to convert it to an ODE system
#This file should be problem specific because we may invert different parameters 

# Problem-specific function to calculate the RHS of ODE: 
# dQdt (= dhdt dq_xdt dq_ydt)
# Q = [h, q_x, q_y] is the solution vector. q_x = h*u_x, q_y = h*u_y
# Q_ghost = [h_ghost, q_x_ghost, q_y_ghost] is the ghost cell values


# In this case, we will invert para = [zb_cells].
function swe_2D_rhs!(dQdt, Q, Q_ghost, para, t, my_mesh_2D, swe_2d_constants, ManningN_cells, ManningN_ghostCells,
    nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
    inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length,
    inletQ_TotalA, inletQ_DryWet,  
    nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
    exitH_internalCellIDs, exitH_faceCentroids, exitH_WSE,
    nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, 
    wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals
    )

    #mesh data
    numOfCells = my_mesh_2D.numOfCells
    numOfFaces = my_mesh_2D.numOfFaces

    cell_areas = my_mesh_2D.cell_areas
    cell_normals = my_mesh_2D.cell_normals
    face_normals = my_mesh_2D.face_normals  
    face_lengths = my_mesh_2D.face_lengths

    faceLeftCellID_Dict = my_mesh_2D.faceLeftCellID_Dict
    faceRightCellID_Dict = my_mesh_2D.faceRightCellID_Dict

    g = swe_2d_constants.g
    h_small = swe_2d_constants.h_small
    RiemannSolver = swe_2d_constants.RiemannSolver

    h = @view Q_cells[:,1]
    q_x = @view Q_cells[:,2]
    q_y = @view Q_cells[:,3]

    #fluxes on faces: 
    fluxes = zeros(eltype(Q), 3, numOfFaces)
    flux = zeros(eltype(Q), 3)

    #zb at faces and slope at cells 
    zb_faces = zeros(eltype(Q), numOfFaces)
    S0 = zeros(eltype(Q), numOfCells)

    #println("time t = ", t)

    #set the parameter values 
    zb_cells = para
    
    #interpolate zb from cell to face 
    cells_to_faces_scalar!(numOfFaces, my_mesh_2D.faceCells_Dict, zb_cells, zb_faces)

    #compute bed slope at cell centers
    compute_scalar_gradients!(numOfCells, cell_areas, cell_normals, face_lengths, 
        my_mesh_2D.cellNodesCount, my_mesh_2D.cellFacesList, my_mesh_2D.cellNeighbors_Dict, 
        zb_cells, S0)

    #bed slope is negative of zb gradient 
    S0 = -S0

    #Process the boundaries: update ghost cells
    process_inlet_q_boundaries(nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
    inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length,
    inletQ_TotalA, inletQ_DryWet, Q, Q_ghost, swe_2D_constants)

    process_exit_h_boundaries(nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
    exitH_internalCellIDs, exitH_faceCentroids, exitH_WSE, Q, Q_ghost, swe_2D_constants)

    process_wall_boundaries(nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, 
    wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals, Q, Q_ghost)


    #loop through all faces
    @inbounds for iFace in 1:my_mesh_2D.numOfFaces
        left_cellID = faceLeftCellID_Dict[iFace]
        right_cellID = faceRightCellID_Dict[iFace]
        face_boundary_ID = my_mesh_2D.faceBoundaryID_Dict[iFace]

        h_face = [h[left_cellID], h[rigth_cellID]]
        q_x_face = [q_x[left_cellID], q_x[right_cellID]]
        q_y_face = [q_y[left_cellID], q_y[right_cellID]]

        if face_boundary_ID !=0                       #boundary face
            if iFace in inletQ_faceIDs
                bcType = "inletQ"
            elseif iFace in exitH_faceIDs
                bcType = "exitH"
            elseif iFace in wall_faceIDs
                bcType = "wall"
            else
                println("Error: wrong face boundary ID")
                exit(-1)  #exit with an error code of -1
            end
        else
            bcType = "internal"
        end 

        if (RiemannSolver == "Roe")
            println("Not implemented yet")
            exit(-1)  #exit with an error code of -1
        elseif (RiemannSolver == "HLL")
            Riemann_2D_hll!(flux, g, h_face, q_x_face, q_y_face, bcType, h_small)
        elseif (RiemannSolver == "HLLC")
            println("Not implemented yet")
            exit(-1)  #exit with an error code of -1
        else
            println("Wrong choice of RiemannSolver")
            exit(-1)  #exit with an error code of -1
        end

        fluxes[1,iFace] = flux[1] 
        fluxes[2,iFace] = flux[2]
        fluxes[3,iFace] = flux[3]
    end 

    #loop through all cells
    @inbounds for iCell in 1:nCells

        #calcuate the RHS of the ODEs

        #fields.dhdt[iCell] = - (flux_east[1]- flux_west[1]) / dx
        dQdt[iCell,1] = - (fluxes[:,iCell+1][1]- fluxes[:,iCell][1]) / dx

        if (h[iCell] <= h_small) #if a dry cell, no flow resistance term
            #fields.dqdt[iCell] = -(flux_east[2] - flux_west[2]) / dx + g * h[iCell] * fields.S0[iCell]
            dQdt[iCell,2] = -(fluxes[:,iCell+1][2] - fluxes[:,iCell][2]) / dx + g * h[iCell] * S0[iCell]
        else
            #fields.dqdt[iCell] = (- (flux_east[2] - flux_west[2]) / dx + g * h[iCell] * fields.S0[iCell]
            #           - g*ManningN^2/max(h[iCell], h_small)^(7.0/3.0)*abs(fields.q[iCell])*fields.q[iCell])
            dQdt[iCell,2] = (- (fluxes[:,iCell+1][2] - fluxes[:,iCell][2]) / dx + g * h[iCell] * S0[iCell]
                       - g*ManningN^2/max(h[iCell], h_small)^(7.0/3.0)*abs(q[iCell])*q[iCell])
        end
    end 

end
