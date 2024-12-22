#semin-discretize the 2D SWE to convert it to an ODE system
#This file should be problem specific because we may invert different parameters 

# Problem-specific function to calculate the RHS of ODE: 
# dQdt (= dhdt dq_xdt dq_ydt)
# Q = [h, q_x, q_y] is the solution vector. q_x = h*u_x, q_y = h*u_y
# Q_ghost = [h_ghost, q_x_ghost, q_y_ghost] is the ghost cell values

using ForwardDiff
using Zygote

# In this case, we will invert para = [zb_cells].
function swe_2d_rhs(Q, Q_ghost, para, t, my_mesh_2D, zb_cells, zb_ghostCells, zb_faces, S0,
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

    #Zygote.ignore() do  
    #    println("within swe_2D_rhs, t =", t)
    #end

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
     #println("zb_ghostCells = ", zb_ghostCells.values)
     #println("S0 = ", S0)
 
     # Process boundaries: update ghost cells values
     # Each boundary treatment function works on different part of Q_ghost. So it is not necessary to create Q_ghost in each of the 
     # boundary treatment function. We can create Q_ghost once and pass it to each of the boundary treatment functions.
     #Q_ghost = zeros(eltype(para), my_mesh_2D.numOfAllBounaryFaces, 3)

     #or write a function to update Q_ghost for all boundaries at once to reduce memory usage
     Q_ghost = process_all_boundaries(Q, my_mesh_2D, nInletQ_BCs, nExitH_BCs, nWall_BCs, nSymm_BCs, inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices,
        inletQ_faceIDs, exitH_faceIDs, wall_faceIDs, symm_faceIDs, inletQ_ghostCellIDs, exitH_ghostCellIDs, wall_ghostCellIDs, symm_ghostCellIDs,
        inletQ_internalCellIDs, exitH_internalCellIDs, wall_internalCellIDs, symm_internalCellIDs, inletQ_faceCentroids, exitH_faceCentroids, wall_faceCentroids, symm_faceCentroids,
        inletQ_faceOutwardNormals, exitH_faceOutwardNormals, wall_outwardNormals, symm_outwardNormals, inletQ_TotalQ, exitH_WSE, wall_H, symm_H,
        inletQ_A, exitH_A, wall_A, symm_A, ManningN_cells, swe_2d_constants)

    #  if nInletQ_BCs > 0
    #     Q_ghost = process_inlet_q_boundaries!(nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs,
    #                                       inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals,
    #                                       inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length,
    #                                       inletQ_TotalA, inletQ_DryWet, Q, Q_ghost, ManningN_cells, swe_2d_constants)
    #  end    

    #  if nExitH_BCs > 0
    #     Q_ghost = process_exit_h_boundaries!(nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs,
    #                                      exitH_internalCellIDs, exitH_faceCentroids, exitH_WSE,
    #                                      Q, Q_ghost, swe_2d_constants)
    #  end

    #  if nWall_BCs > 0
    #     Q_ghost = process_wall_boundaries!(nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs,
    #                                       wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals,
    #                                       Q, Q_ghost)
    #  end

    #  if nSymm_BCs > 0
    #     Q_ghost = process_symmetry_boundaries!(nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs,
    #                                       symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals,
    #                                       Q, Q_ghost)
    #  end
     
     #println("within swe_2D_rhs!")
     #println("Q_cells value = ", ForwardDiff.value.(Q))
     #println("Q_cells partials = ", ForwardDiff.partials.(Q))
     #println("Q_ghost value = ", ForwardDiff.value.(Q_ghost))
     #println("Q_ghost partials = ", ForwardDiff.partials.(Q_ghost))
     
     #throw("stop here")

     # Initialize fluxes
     #dQdt = zeros(eltype(Q), size(Q))
 
     # Loop through all cells to calculate the fluxes on faces
     updates = map(1:numOfCells) do iCell           # .= is in-place mutation!
         cell_area = cell_areas[iCell]
 
         # Initialize flux accumulation
         flux_sum = zeros(eltype(Q), 3)
         
         for iFace in 1:my_mesh_2D.cellNodesCount[iCell]
             faceID = abs(my_mesh_2D.cellFacesList[iCell,:][iFace])
             left_cellID = iCell
             right_cellID = abs(cellNeighbors_Dict[iCell][iFace])
             
             faceBoundaryID = my_mesh_2D.faceBoundaryID_Dict[faceID]
             face_normal = cell_normals[iCell][iFace]
             
             if faceBoundaryID == 0  # internal face
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
             flux_sum = flux_sum .+ flux .* face_lengths[faceID]
         end

         if iCell == -1
            println("flux_sum value = ", ForwardDiff.value.(flux_sum))
            println("flux_sum partials = ", ForwardDiff.partials.(flux_sum))
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
             #u_mag = sqrt(u_temp^2 + v_temp^2)
             u_mag = max(sqrt(u_temp^2 + v_temp^2), sqrt(eps(eltype(u_temp))))
             
             friction_x = g * ManningN_cells[iCell]^2 / h[iCell]^(1.0/3.0) * u_mag * u_temp
             friction_y = g * ManningN_cells[iCell]^2 / h[iCell]^(1.0/3.0) * u_mag * v_temp

             if iCell == -1
                println("u_temp = ", u_temp)
                println("v_temp = ", v_temp)
                println("u_mag = ", u_mag)
                println("ManningN_cells[iCell] = ", ManningN_cells[iCell])
                println("q_x[iCell] = ", q_x[iCell])
                println("q_y[iCell] = ", q_y[iCell])
                println("h[iCell] = ", h[iCell])
                println("S0[iCell, 1] = ", S0[iCell, 1])
                println("S0[iCell, 2] = ", S0[iCell, 2])
                println("cell_area = ", cell_area)
                println("friction_x value = ", ForwardDiff.value.(friction_x))
                println("friction_x partials = ", ForwardDiff.partials.(friction_x))
                println("friction_y value = ", ForwardDiff.value.(friction_y))
                println("friction_y partials = ", ForwardDiff.partials.(friction_y))
             end
 
             source_terms .= [zero(eltype(Q)), 
                              (g * h[iCell] * S0[iCell, 1] - friction_x) * cell_area,
                              (g * h[iCell] * S0[iCell, 2] - friction_y) * cell_area]
         end

         if iCell == -1
            println("source_terms value = ", ForwardDiff.value.(source_terms))
            println("source_terms partials = ", ForwardDiff.partials.(source_terms))
         end
 
         # Return the update for this cell (without in-place mutation)
         (-flux_sum .+ source_terms) ./ cell_area
     end

     #convert updates to a 2D array: vcat organizes vectors as rows
     # Stacks vectors vertically, treating each vector as a row of the resulting matrix.
     # The transpose (') ensures the vectors are treated as rows when concatenated.
     dQdt = vcat(updates'...)

     #println("dQdt value = ", ForwardDiff.value.(dQdt))
     #println("dQdt partials = ", ForwardDiff.partials.(dQdt))

     #println("dQ/dpara = ", ForwardDiff.partials(Q))
     #throw("stop here")
 
     return dQdt
end
