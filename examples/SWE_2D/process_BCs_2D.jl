#preprcess all boundary conditions, such as calculate the list of faces, internal cells, ghost cells, etc.
#These functions are supposed to be called only once before the time loop

using LinearAlgebra
using ForwardDiff
using Zygote




#process all boundaries in just one function call to update Q_ghostCells: called every time step to update the boundary condition
#return a new Q_ghostCells
function process_all_boundaries_2d(Q_cells, my_mesh_2D, nInletQ_BCs, nExitH_BCs, nWall_BCs, nSymm_BCs, inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices,
    inletQ_faceIDs, exitH_faceIDs, wall_faceIDs, symm_faceIDs, inletQ_ghostCellIDs, exitH_ghostCellIDs, wall_ghostCellIDs, symm_ghostCellIDs,
    inletQ_internalCellIDs, exitH_internalCellIDs, wall_internalCellIDs, symm_internalCellIDs, inletQ_faceCentroids, exitH_faceCentroids, wall_faceCentroids, symm_faceCentroids,
    inletQ_faceOutwardNormals, exitH_faceOutwardNormals, wall_outwardNormals, symm_outwardNormals, inletQ_TotalQ, exitH_WSE, wall_H, symm_H,
    inletQ_A, exitH_A, wall_A, symm_A, ManningN_cells, swe_2D_constants)

    h = @view Q_cells[:, 1]         #view just creates a reference to the data, not a copy
    q_x = @view Q_cells[:, 2]
    q_y = @view Q_cells[:, 3]

    h_ghost = zeros(eltype(para), my_mesh_2D.numOfAllBounaryFaces)
    q_x_ghost = zeros(eltype(para), my_mesh_2D.numOfAllBounaryFaces)
    q_y_ghost = zeros(eltype(para), my_mesh_2D.numOfAllBounaryFaces)

    #h_ghost = copy(Q_ghostCells[:, 1])
    #q_x_ghost = copy(Q_ghostCells[:, 2])
    #q_y_ghost = copy(Q_ghostCells[:, 3])

    #for inlet-q boundaries
    # Loop through all inlet-q boundaries
    for iInletQ in 1:nInletQ_BCs
        current_boundaryFaceIDs = inletQ_faceIDs[iInletQ]
        current_ghostCellIDs = inletQ_ghostCellIDs[iInletQ]
        current_internalCellIDs = inletQ_internalCellIDs[iInletQ]
        current_inletQ_Length = inletQ_Length[iInletQ]
        current_inletQ_TotalQ = inletQ_TotalQ[iInletQ]

        # Create new arrays for this boundary's updates
        h_new = zeros(eltype(Q_cells), length(current_ghostCellIDs))
        q_x_new = zeros(eltype(Q_cells), length(current_ghostCellIDs))
        q_y_new = zeros(eltype(Q_cells), length(current_ghostCellIDs))

        # First pass: compute total_A and dry/wet status
        total_A = zero(eltype(Q_cells))
        current_inletQ_DryWet = zeros(Int, length(current_ghostCellIDs))

        for (iFace, internalCellID) in enumerate(current_internalCellIDs)
            ManningN_face = ManningN_cells[internalCellID]

            if h[internalCellID] > swe_2D_constants.h_small
                current_inletQ_DryWet[iFace] = 1
                total_A = total_A + current_inletQ_Length[iFace]^(5.0 / 3.0) * h[internalCellID] / ManningN_face
            end
        end

        if total_A <= 1e-10
            error("Total cross-sectional conveyance for inlet-q boundary $iInletQ is not positive: $total_A")
        end

        # Second pass: compute ghost cell values
        for (iFace, (ghostCellID, internalCellID)) in enumerate(zip(current_ghostCellIDs, current_internalCellIDs))
            ManningN_face = ManningN_cells[internalCellID]
            h_new[iFace] = h[internalCellID]

            if current_inletQ_DryWet[iFace] == 0
                q_x_new[iFace] = zero(eltype(Q_cells))
                q_y_new[iFace] = zero(eltype(Q_cells))
            else
                face_normal = inletQ_faceOutwardNormals[iInletQ][iFace, :]
                velocity_normal = (current_inletQ_TotalQ / total_A *
                                   current_inletQ_Length[iFace]^(2.0 / 3.0) / ManningN_face)

                q_x_new[iFace] = -h_new[iFace] * velocity_normal * face_normal[1]
                q_y_new[iFace] = -h_new[iFace] * velocity_normal * face_normal[2]
            end
        end

        # Update ghost arrays using update_1d_array
        h_ghost = update_1d_array(h_ghost, current_ghostCellIDs, h_new)
        q_x_ghost = update_1d_array(q_x_ghost, current_ghostCellIDs, q_x_new)
        q_y_ghost = update_1d_array(q_y_ghost, current_ghostCellIDs, q_y_new)
    end

    #for exit-h boundaries
    #loop through all exit-h boundaries
    for iExitH in 1:nExitH_BCs
        #iBoundary = exitH_BC_indices[iExitH]

        current_ghostCellIDs = exitH_ghostCellIDs[iExitH]  # ghost cell IDs for this boundary
        current_internalCellIDs = exitH_internalCellIDs[iExitH]  # internal cell IDs for this boundary
        current_faceCentroids = exitH_faceCentroids[iExitH]  # face centroids for this boundary

        # Calculate new h_ghost values
        h_new = convert.(eltype(Q_cells), max.(swe_2D_constants.h_small,
            exitH_WSE[iExitH] .- current_faceCentroids[:, 3]))

        # Update arrays using update_1d_array
        h_ghost = update_1d_array(h_ghost, current_ghostCellIDs, h_new)
        q_x_ghost = update_1d_array(q_x_ghost, current_ghostCellIDs, q_x[current_internalCellIDs])
        q_y_ghost = update_1d_array(q_y_ghost, current_ghostCellIDs, q_y[current_internalCellIDs])
    end

    #for wall boundaries
    #loop through all wall boundaries
    for iWall in 1:nWall_BCs
        # Get indices
        current_ghostCellIDs = wall_ghostCellIDs[iWall]
        current_internalCellIDs = wall_internalCellIDs[iWall]

        # updates (without in-place mutation)
        h_ghost = update_1d_array(h_ghost, current_ghostCellIDs, h[current_internalCellIDs])
        q_x_ghost = update_1d_array(q_x_ghost, current_ghostCellIDs, -q_x[current_internalCellIDs])
        q_y_ghost = update_1d_array(q_y_ghost, current_ghostCellIDs, -q_y[current_internalCellIDs])
    end

    #for symmetry boundaries
    #loop through all symmetry boundaries
    for iSymm in 1:nSymm_BCs
        current_ghostCellIDs = symm_ghostCellIDs[iSymm]
        current_internalCellIDs = symm_internalCellIDs[iSymm]
        face_normals = symm_outwardNormals[iSymm]

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

