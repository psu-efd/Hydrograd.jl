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

     # Mesh data
     numOfCells = my_mesh_2D.numOfCells
     maxNumOfCellFaces = my_mesh_2D.maxNumOfCellFaces
     
     cell_areas = my_mesh_2D.cell_areas
     cell_normals = my_mesh_2D.cell_normals
     face_lengths = my_mesh_2D.face_lengths
     cellNeighbors_Dict = my_mesh_2D.cellNeighbors_Dict
     
     g = swe_2d_constants.g
     h_small = swe_2d_constants.h_small
     RiemannSolver = swe_2d_constants.RiemannSolver
     
     h = @view Q[:, 1]  # water depth
     q_x = @view Q[:, 2]  # discharge in x
     q_y = @view Q[:, 3]  # discharge in y
 
     # Set parameter values
     zb_cells_current = para  # para.clone()
     
     # Interpolate zb from cell to face and compute bed slope at cells
     zb_ghostCells, zb_faces, S0 = interploate_zb_from_cell_to_face_and_compute_S0(my_mesh_2D, zb_cells_current)
 
     # Process boundaries: update ghost cells
     # Each boundary treatment function works on different part of Q_ghost. So it is not necessary to create Q_ghost in each of the 
     # boundary treatment function. We can create Q_ghost once and pass it to each of the boundary treatment functions.
     #Q_ghost = zeros(eltype(Q), size(Q))
     Q_ghost = process_inlet_q_boundaries(nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs,
                                          inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals,
                                          inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length,
                                          inletQ_TotalA, inletQ_DryWet, Q, Q_ghost, ManningN_cells, swe_2d_constants)
     
     Q_ghost = process_exit_h_boundaries(nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs,
                                         exitH_internalCellIDs, exitH_faceCentroids, exitH_WSE,
                                         Q, Q_ghost, swe_2d_constants)
     
     Q_ghost = process_wall_boundaries(nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs,
                                       wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals,
                                       Q, Q_ghost)
     
     Q_ghost = process_symmetry_boundaries(nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs,
                                           symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals,
                                           Q, Q_ghost)
     
     # Initialize fluxes
     #dQdt = zeros(eltype(Q), size(Q))
 
     # Loop through all cells to calculate the fluxes on faces
     for iCell in 1:numOfCells
         cell_area = cell_areas[iCell]
 
         # Initialize flux accumulation
         flux_sum = zeros(eltype(Q), 3)
         
         for iFace in 1:my_mesh_2D.cellNodesCount[iCell]
             faceID = abs(my_mesh_2D.cellFacesList[iCell,:][iFace])
             left_cellID = iCell
             right_cellID = abs(cellNeighbors_Dict[iCell][iFace])
             
             faceBoundaryID = my_mesh_2D.faceBoundaryID_Dict[faceID]
             face_normal = cell_normals[iCell][iFace]
             
             if faceBoundaryID == -1  # internal face
                 hL, huL, hvL = h[left_cellID], q_x[left_cellID], q_y[left_cellID]
                 hR, huR, hvR = h[right_cellID], q_x[right_cellID], q_y[right_cellID]
             else  # boundary face
                 hL, huL, hvL = h[left_cellID], q_x[left_cellID], q_y[left_cellID]
                 hR = Q_ghost[right_cellID, 1]
                 huR = Q_ghost[right_cellID, 2]
                 hvR = Q_ghost[right_cellID, 3]
             end
             
             if RiemannSolver == "Roe"
                 flux = Riemann_2D_Roe(hL, huL, hvL, hR, huR, hvR, g, face_normal, hmin=h_small)
             elseif RiemannSolver == "HLL"
                 error("HLL solver not implemented yet")
             elseif RiemannSolver == "HLLC"
                 error("HLLC solver not implemented yet")
             else
                 error("Wrong choice of RiemannSolver")
             end
             
             # Accumulate flux contribution
             flux_sum .= flux_sum .+ flux .* face_lengths[faceID]
         end
 
         # Source terms
         source_terms = zeros(eltype(Q), 3)
         if h[iCell] <= h_small
             source_terms .= [eltype(Q)(0.0), 
                              g * h[iCell] * S0[iCell, 1] * cell_area, 
                              g * h[iCell] * S0[iCell, 2] * cell_area]
         else
             u_temp = q_x[iCell] / h[iCell]
             v_temp = q_y[iCell] / h[iCell]
             u_mag = sqrt(u_temp^2 + v_temp^2)
             
             friction_x = g * ManningN_cells[iCell]^2 / h[iCell]^(1/3) * u_mag * u_temp
             friction_y = g * ManningN_cells[iCell]^2 / h[iCell]^(1/3) * u_mag * v_temp
 
             source_terms .= [0.0, 
                              (g * h[iCell] * S0[iCell, 1] - friction_x) * cell_area,
                              (g * h[iCell] * S0[iCell, 2] - friction_y) * cell_area]
         end
 
         # Compute final contribution to dQdt (combine flux and source terms)
         dQdt[iCell, :] .= (-flux_sum .+ source_terms) ./ cell_area
     end
 
     return dQdt
end
