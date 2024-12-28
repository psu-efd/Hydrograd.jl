#preprcess all boundary conditions, such as calculate the list of faces, internal cells, ghost cells, etc.
#These functions are supposed to be called only once before the time loop

using LinearAlgebra
using ForwardDiff
using Zygote

#process all boundaries in just one function call to update Q_ghostCells: called every time step to update the boundary condition
#return a new Q_ghostCells
function process_all_boundaries_2d(settings, Q_cells, my_mesh_2D, boundary_conditions, ManningN_cells, zb_faces, swe_2D_constants, inletQ_Length, inletQ_TotalQ, exitH_WSE)

    h = @view Q_cells[:, 1]         
    q_x = @view Q_cells[:, 2]
    q_y = @view Q_cells[:, 3]

    h_ghost = zeros(eltype(Q_cells), my_mesh_2D.numOfAllBounaryFaces)
    q_x_ghost = zeros(eltype(Q_cells), my_mesh_2D.numOfAllBounaryFaces)
    q_y_ghost = zeros(eltype(Q_cells), my_mesh_2D.numOfAllBounaryFaces)

    #for inlet-q boundaries
    # Loop through all inlet-q boundaries
    for iInletQ in 1:boundary_conditions.nInletQ_BCs

        iBoundary = boundary_conditions.inletQ_BC_indices[iInletQ]

        if settings.bVerbose
            Zygote.ignore() do
                println("Processing INLET-Q boundary ", iInletQ, " with index in BC list ", iBoundary)
            end
        end

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
        #current_inletQ_DryWet = map(internalCellID -> h[internalCellID] > swe_2D_constants.h_small ? 1 : 0, current_internalCellIDs)  #does this break AD?
        current_inletQ_DryWet = map(internalCellID -> 1, current_internalCellIDs)

        # Compute total_A only for wet faces
        #total_A = sum(current_inletQ_Length[iFace]^(5.0 / 3.0) * h[current_internalCellIDs[iFace]] / ManningN_cells[current_internalCellIDs[iFace]]
        #              for iFace in 1:length(current_internalCellIDs) if current_inletQ_DryWet[iFace] == 1
        #)

        #contributions = [
        #    current_inletQ_Length[iFace]^(5.0 / 3.0) * h[current_internalCellIDs[iFace]] / ManningN_cells[current_internalCellIDs[iFace]]
        #    for iFace in 1:length(current_internalCellIDs) if current_inletQ_DryWet[iFace] == 1
        #]

        # contributions = map(1:length(current_internalCellIDs)) do iFace
        #     base_contribution = current_inletQ_Length[iFace]^(5.0 / 3.0) * 
        #                       h[current_internalCellIDs[iFace]] / 
        #                       ManningN_cells[current_internalCellIDs[iFace]]
        #     # Multiply by wet/dry status (0 or 1)
        #     base_contribution * Float64(current_inletQ_DryWet[iFace])
        # end

        #fake contributions for debugging
        contributions = [33.33333333322568, 33.333333333285225, 33.333333333489115]

        # Debug type information
        Zygote.ignore() do
            @show eltype(contributions)
            @show typeof(contributions)
            all(x -> isa(x, Float64), contributions) || 
                @warn "Non-Float64 elements found in contributions"

            @show size(contributions)
            @show contributions
        end

        total_A = sum(contributions)

        @assert total_A > 1e-10 "Total cross-sectional conveyance for inlet-q boundary $iInletQ is not positive: $total_A"

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

        # Calculate q_updates in a non-mutating way
        # q_x_new = map(enumerate(current_internalCellIDs)) do (iFace, internalCellID)
        #     ManningN_face = ManningN_cells[internalCellID]
        #     face_normal = boundary_conditions.inletQ_faceOutwardNormals[iInletQ][iFace, :]
        #     velocity_normal = current_inletQ_TotalQ / total_A *
        #                     current_inletQ_Length[iFace]^(2.0 / 3.0) / ManningN_face
            
        #     # Use multiplication by wet/dry status instead of if/else
        #     factor = Float64(current_inletQ_DryWet[iFace])
        #     -h[internalCellID] * velocity_normal * face_normal[1] * factor
        # end

        Zygote.ignore() do
            @show typeof(current_inletQ_TotalQ)
            @show typeof(total_A)
            @show typeof(current_inletQ_Length)
            @show typeof(ManningN_cells)

            @show current_inletQ_TotalQ
            @show total_A
            @show current_inletQ_Length
            @show ManningN_cells
        end

        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (current_inletQ_TotalQ / total_A * current_inletQ_Length[iFace]^(2.0 / 3.0) / ManningN_cells[internalCellID]) 
        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (current_inletQ_TotalQ / total_A / ManningN_cells[internalCellID]) #break
        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (current_inletQ_TotalQ / ManningN_cells[internalCellID]) #break
        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (current_inletQ_TotalQ / total_A ) #works
        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * ( ManningN_cells[internalCellID]) #breaks
        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (current_inletQ_Length[iFace]^(2.0 / 3.0)) #works; current_inletQ_Length from a static struct (?)
        #                         * boundary_conditions.inletQ_faceOutwardNormals[iInletQ][iFace, 1] * Float64(current_inletQ_DryWet[iFace]), enumerate(current_internalCellIDs))
        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (current_inletQ_TotalQ  * current_inletQ_Length[iFace]^(2.0 / 3.0) / ManningN_cells[internalCellID])  #break
        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (current_inletQ_TotalQ  * current_inletQ_Length[iFace]^(2.0 / 3.0) / ManningN_cells[internalCellID])  #break
        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (current_inletQ_TotalQ )    #works
        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (current_inletQ_TotalQ /2.0 )    #works
        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (total_A )    #break
        #q_x_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (current_inletQ_TotalQ / ManningN_cells[internalCellID])  #break
        #                 * boundary_conditions.inletQ_faceOutwardNormals[iInletQ][iFace, 1] * Float64(current_inletQ_DryWet[iFace]), enumerate(current_internalCellIDs))

        # q_y_new = map(enumerate(current_internalCellIDs)) do (iFace, internalCellID)
        #     ManningN_face = ManningN_cells[internalCellID]
        #     face_normal = boundary_conditions.inletQ_faceOutwardNormals[iInletQ][iFace, :]
        #     velocity_normal = current_inletQ_TotalQ / total_A *
        #                     current_inletQ_Length[iFace]^(2.0 / 3.0) / ManningN_face
            
        #     # Use multiplication by wet/dry status instead of if/else
        #     factor = Float64(current_inletQ_DryWet[iFace])
        #     -h[internalCellID] * velocity_normal * face_normal[2] * factor
        # end

        #q_y_new = map(((iFace, internalCellID),) -> -h[internalCellID] * (current_inletQ_TotalQ / total_A * current_inletQ_Length[iFace]^(2.0 / 3.0) / ManningN_cells[internalCellID]) 
        #                 * boundary_conditions.inletQ_faceOutwardNormals[iInletQ][iFace, 2] * Float64(current_inletQ_DryWet[iFace]), enumerate(current_internalCellIDs))
        #q_y_new = map(((iFace, internalCellID),) -> -h[internalCellID] * 0.0 
        #                         * boundary_conditions.inletQ_faceOutwardNormals[iInletQ][iFace, 2] * Float64(current_inletQ_DryWet[iFace]), enumerate(current_internalCellIDs))

        Zygote.ignore() do
            @show typeof(q_x_new)
            @show typeof(q_y_new)
            @show size(q_x_new)
            @show size(q_y_new)
            @show q_x_new
            @show q_y_new
        end

        #q_x_new = convert.(eltype(Q_cells), q_x_new)
        #q_y_new = convert.(eltype(Q_cells), q_y_new)

        #fake q_x_new and q_y_new for debugging
        #q_x_new = [0.0 for iFace in 1:length(current_internalCellIDs)]
        #q_y_new = [0.0 for iFace in 1:length(current_internalCellIDs)]

        # Update ghost arrays using update_1d_array
        h_ghost = update_1d_array(h_ghost, current_ghostCellIDs, h_new)
        q_x_ghost = update_1d_array(q_x_ghost, current_ghostCellIDs, q_x_new)
        q_y_ghost = update_1d_array(q_y_ghost, current_ghostCellIDs, q_y_new)
    end

    #for exit-h boundaries
    #loop through all exit-h boundaries
    for iExitH in 1:boundary_conditions.nExitH_BCs
        iBoundary = boundary_conditions.exitH_BC_indices[iExitH]

        if settings.bVerbose
            Zygote.ignore() do
                println("Processing EXIT-H boundary ", iExitH, " with index in BC list ", iBoundary)
            end
        end

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
        iBoundary = boundary_conditions.wall_BC_indices[iWall]

        if settings.bVerbose
            Zygote.ignore() do
                println("Processing WALL boundary ", iWall, " with index in BC list ", iBoundary)
            end
        end

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
        iBoundary = boundary_conditions.symm_BC_indices[iSymm]

        if settings.bVerbose
            Zygote.ignore() do
                println("Processing SYMMETRY boundary ", iSymm, " with index in BC list ", iBoundary)
            end
        end

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

    Zygote.ignore() do
        @show size(new_Q_ghostCells)
        @show typeof(new_Q_ghostCells)
        @show new_Q_ghostCells
    end

    return new_Q_ghostCells
end

#update inletQ_TotalQ
function update_inletQ_TotalQ(inlet_discharges_current)
    #update inletQ_TotalQ based on the provided inlet_discharges_current

    inletQ_TotalQ =  [inlet_discharges_current[iInletQ] for iInletQ in 1:(length(inlet_discharges_current))]

    return inletQ_TotalQ

end

