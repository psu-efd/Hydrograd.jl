#preprcess all boundary conditions, such as calculate the list of faces, internal cells, ghost cells, etc.
#These functions are supposed to be called only once before the time loop

using LinearAlgebra
using ForwardDiff
using Zygote

#process all boundaries in just one function call to update Q_ghostCells: called every time step to update the boundary condition
#return a new Q_ghostCells
function process_all_boundaries_2d(Q_cells, my_mesh_2D, boundary_conditions, ManningN_cells, zb_faces, swe_2D_constants, inletQ_Length, inletQ_TotalQ, exitH_WSE)

    h = copy(Q_cells[:, 1])         
    q_x = copy(Q_cells[:, 2])
    q_y = copy(Q_cells[:, 3])

    h_ghost = zeros(eltype(para), my_mesh_2D.numOfAllBounaryFaces)
    q_x_ghost = zeros(eltype(para), my_mesh_2D.numOfAllBounaryFaces)
    q_y_ghost = zeros(eltype(para), my_mesh_2D.numOfAllBounaryFaces)

    #for inlet-q boundaries
    # Loop through all inlet-q boundaries
    for iInletQ in 1:boundary_conditions.nInletQ_BCs
        current_boundaryFaceIDs = boundary_conditions.inletQ_faceIDs[iInletQ]
        current_ghostCellIDs = boundary_conditions.inletQ_ghostCellIDs[iInletQ]
        current_internalCellIDs = boundary_conditions.inletQ_internalCellIDs[iInletQ]

        current_inletQ_Length = inletQ_Length[iInletQ]
        current_inletQ_TotalQ = inletQ_TotalQ[iInletQ]

        # Create new arrays for this boundary's updates
        h_new = zeros(eltype(Q_cells), length(current_ghostCellIDs))
        q_x_new = zeros(eltype(Q_cells), length(current_ghostCellIDs))
        q_y_new = zeros(eltype(Q_cells), length(current_ghostCellIDs))

        # First pass: compute total_A and dry/wet status
        # current_inletQ_DryWet: 1 for wet, 0 for dry
        current_inletQ_DryWet = map(internalCellID -> h[internalCellID] > swe_2D_constants.h_small ? 1 : 0, current_internalCellIDs)

        # Compute total_A only for wet faces
        total_A = sum(current_inletQ_Length[iFace]^(5.0 / 3.0) * h[current_internalCellIDs[iFace]] / ManningN_cells[current_internalCellIDs[iFace]]
                      for iFace in 1:length(current_internalCellIDs) if current_inletQ_DryWet[iFace] == 1
        )

        if total_A <= 1e-10
            error("Total cross-sectional conveyance for inlet-q boundary $iInletQ is not positive: $total_A")
        end

        # Second pass: compute ghost cell values
        h_new = map(internalCellID -> h[internalCellID], current_internalCellIDs)

        q_updates = [(
           let
               ManningN_face = ManningN_cells[internalCellID]
               face_normal = boundary_conditions.inletQ_faceOutwardNormals[iInletQ][iFace, :]
               velocity_normal = (current_inletQ_TotalQ / total_A *
                                current_inletQ_Length[iFace]^(2.0 / 3.0) / ManningN_face)
                if current_inletQ_DryWet[iFace] == 0
                   (eltype(Q_cells)(0.0), eltype(Q_cells)(0.0))
               else
                   (-h[internalCellID] * velocity_normal * face_normal[1],
                    -h[internalCellID] * velocity_normal * face_normal[2])
               end
           end
       ) for (iFace, internalCellID) in enumerate(current_internalCellIDs)]

        q_x_new = [q_updates[iFace][1] for iFace in 1:length(current_internalCellIDs)]
        q_y_new = [q_updates[iFace][2] for iFace in 1:length(current_internalCellIDs)]

        # Update ghost arrays using update_1d_array
        h_ghost = update_1d_array(h_ghost, current_ghostCellIDs, h_new)
        q_x_ghost = update_1d_array(q_x_ghost, current_ghostCellIDs, q_x_new)
        q_y_ghost = update_1d_array(q_y_ghost, current_ghostCellIDs, q_y_new)
    end

    #for exit-h boundaries
    #loop through all exit-h boundaries
    for iExitH in 1:boundary_conditions.nExitH_BCs
        #iBoundary = exitH_BC_indices[iExitH]

        current_ghostCellIDs = boundary_conditions.exitH_ghostCellIDs[iExitH]  # ghost cell IDs for this boundary
        current_internalCellIDs = boundary_conditions.exitH_internalCellIDs[iExitH]  # internal cell IDs for this boundary
        current_faceIDs = boundary_conditions.exitH_faceIDs[iExitH]  # face IDs for this boundary

        #current_faceCentroids = exitH_faceCentroids[iExitH]  # face centroids for this boundary
        
        # Calculate new h_ghost values
        h_new = convert.(eltype(Q_cells), max.(swe_2D_constants.h_small,
            exitH_WSE[iExitH] .- zb_faces[current_faceIDs]))

        # Update arrays using update_1d_array
        h_ghost = update_1d_array(h_ghost, current_ghostCellIDs, h_new)
        q_x_ghost = update_1d_array(q_x_ghost, current_ghostCellIDs, q_x[current_internalCellIDs])
        q_y_ghost = update_1d_array(q_y_ghost, current_ghostCellIDs, q_y[current_internalCellIDs])
    end

    #for wall boundaries
    #loop through all wall boundaries
    for iWall in 1:boundary_conditions.nWall_BCs
        # Get indices
        current_ghostCellIDs = boundary_conditions.wall_ghostCellIDs[iWall]
        current_internalCellIDs = boundary_conditions.wall_internalCellIDs[iWall]

        # updates (without in-place mutation)
        h_ghost = update_1d_array(h_ghost, current_ghostCellIDs, h[current_internalCellIDs])
        q_x_ghost = update_1d_array(q_x_ghost, current_ghostCellIDs, -q_x[current_internalCellIDs])
        q_y_ghost = update_1d_array(q_y_ghost, current_ghostCellIDs, -q_y[current_internalCellIDs])
    end

    #for symmetry boundaries
    #loop through all symmetry boundaries
    for iSymm in 1:boundary_conditions.nSymm_BCs
        current_ghostCellIDs = boundary_conditions.symm_ghostCellIDs[iSymm]
        current_internalCellIDs = boundary_conditions.symm_internalCellIDs[iSymm]
        face_normals = boundary_conditions.symm_outwardNormals[iSymm]

        # Compute v·n (dot product of velocity and face normal)
        v_dot_n = q_x[current_internalCellIDs] .* face_normals[:, 1] .+
                  q_y[current_internalCellIDs] .* face_normals[:, 2]

        # Update depths using update_1d_array
        h_ghost = update_1d_array(h_ghost, current_ghostCellIDs, h[current_internalCellIDs])

        # Update velocities using update_1d_array
        q_x_new = q_x[current_internalCellIDs] .- 2.0 .* v_dot_n .* face_normals[:, 1]
        q_y_new = q_y[current_internalCellIDs] .- 2.0 .* v_dot_n .* face_normals[:, 2]

        q_x_ghost = update_1d_array(q_x_ghost, current_ghostCellIDs, q_x_new)
        q_y_ghost = update_1d_array(q_y_ghost, current_ghostCellIDs, q_y_new)
    end

    # Stack the updated ghost cell values and return new Q_ghostCells
    new_Q_ghostCells = hcat(h_ghost, q_x_ghost, q_y_ghost)
    return new_Q_ghostCells
end

