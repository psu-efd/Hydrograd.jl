#semin-discretize the 2D SWE to convert it to an ODE system
#This file should be problem specific because we may invert different parameters 

# Problem-specific function to calculate the RHS of ODE: 
# dQdt (= dhdt dq_xdt dq_ydt)
# Q = [h, q_x, q_y] is the solution vector. q_x = h*u_x, q_y = h*u_y
# Q_ghost = [h_ghost, q_x_ghost, q_y_ghost] is the ghost cell values

# Arguments with "_passed" are passed to the function. They are only 
# used for forward simulations. For inversion and sensitivity analysis, 
# ManningN_cells, ManningN_ghostCells, inletQ_TotalQ, exitH_WSE, 
# zb_cells, zb_ghostCells, zb_faces, and S0 are all derived from params_array.

using ForwardDiff
using Zygote


#function swe_2d_rhs(Q, params_array, active_range, param_ranges, t, settings,
#    my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants, inletQ_Length, exitH_WSE,
#    )

function swe_2d_rhs(Q, params_array, active_range, param_ranges, t, settings,
    my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants,
    ManningN_cells_passed, ManningN_ghostCells_passed, inletQ_Length, inletQ_TotalQ_passed, exitH_WSE_passed,
    zb_cells_passed, zb_ghostCells_passed, zb_faces_passed, S0_passed)

    # These are the arguments that are passed to the function. They are only 
    # used for forward simulations. For inversion and sensitivity analysis, 
    # ManningN_cells, ManningN_ghostCells, inletQ_TotalQ, exitH_WSE, 
    # zb_cells, zb_ghostCells, zb_faces, and S0 are all derived from params_array (computed later in the function)
    ManningN_cells_local = ManningN_cells_passed             #create a reference to ManningN_cells_passed. This makes Zygote happy (Don't know why)
    ManningN_ghostCells_local = ManningN_ghostCells_passed
    inletQ_TotalQ_local = inletQ_TotalQ_passed
    exitH_WSE_local = exitH_WSE_passed
    zb_cells_local = zb_cells_passed                #not used for anything
    zb_ghostCells_local = zb_ghostCells_passed
    zb_faces_local = zb_faces_passed
    S0_local = S0_passed

    #get the data type of Q
    data_type = eltype(Q)

    Zygote.ignore() do
        if settings.bVerbose
            #println("within swe_2D_rhs, t =", t)
            #println("asserting data_type = ", data_type)
        end

        @assert data_type <: Real "data_type must be a subtype of Real for AD compatibility"
    end

    # Extract parameters from the 1D array
    zb_cells_current = @view params_array[param_ranges.zb_start:param_ranges.zb_end]
    ManningN_list_current = @view params_array[param_ranges.manning_start:param_ranges.manning_end]

    nInletQ_BCs = srh_all_Dict["nInletQ_BCs"]

    if nInletQ_BCs > 0
        inlet_discharges_current = @view params_array[param_ranges.inletQ_start:param_ranges.inletQ_end]
    else
        inlet_discharges_current = nothing
    end

    Zygote.ignore() do
        if settings.bVerbose

            #@show typeof(Q)
            #@show typeof(params_array)

            #@show typeof(zb_cells_current)
            #@show typeof(ManningN_list_current)
            #@show typeof(inlet_discharges_current)

            #@show zb_cells_current
            #@show ManningN_list_current
            #@show inlet_discharges_current

        end

    end

    #return zeros(size(Q)) #zeros(data_type, size(Q))
    #return Q .* 1.0

    g = swe_2D_constants.g
    h_small = swe_2D_constants.h_small
    RiemannSolver = swe_2D_constants.RiemannSolver

    h = @view Q[:, 1]  # water depth      #view of 2D array does not create a new array; only a reference
    q_x = @view Q[:, 2]  # discharge in x
    q_y = @view Q[:, 3]  # discharge in y
    #h = Q[:, 1]  # water depth        #slice of 2D array creates a new array (no need; waste of memory)
    #q_x = Q[:, 2]  # discharge in x
    #q_y = Q[:, 3]  # discharge in y

    # Make the return value depend on parameters
    #return Q * 1.0 .+ reshape(zb_cells_current, :, 1) .* ones(1, 3)

    # For the case of inversion or sensitivity analysis, and if zb is an active parameter, 
    # we need to interpolate zb from cell to face and compute bed slope at cells
    if (settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis) &&
       "zb" in settings.inversion_settings.active_param_names

        Zygote.ignore() do
            if settings.bVerbose
                #println("calling interploate_zb_from_cell_to_face_and_compute_S0")
            end
        end

        zb_ghostCells_local, zb_faces_local, S0_local = interploate_zb_from_cell_to_face_and_compute_S0(my_mesh_2D, zb_cells_current)
    end

    Zygote.ignore() do
        if settings.bVerbose
            #@show typeof(zb_ghostCells_local)
            #@show typeof(zb_faces_local)
            #@show typeof(S0_local)
            #@show size(S0_local)
            #@show S0_local
        end
    end

    #For the case of inversion or sensitivity analysis, and if ManningN is an active parameter, 
    # we need to update ManningN at cells and ghost cells
    if (settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis) &&
       "ManningN" in settings.inversion_settings.active_param_names

        Zygote.ignore() do
            if settings.bVerbose
                #println("calling update_ManningN")
            end
        end

        ManningN_cells_local, ManningN_ghostCells_local = update_ManningN(my_mesh_2D, srh_all_Dict, ManningN_list_current)
    end

    Zygote.ignore() do
        if settings.bVerbose
            #@show typeof(ManningN_cells_local)
            #@show typeof(ManningN_ghostCells_local)
            #@show size(ManningN_cells_local)
            #@show ManningN_cells_local
        end
    end

    #For the case of inversion or sensitivity analysis, and if Q is an active parameter, 
    # we need to update inletQ_TotalQ based on the provided inlet_discharges_current
    if (settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis) &&
       "Q" in settings.inversion_settings.active_param_names

        Zygote.ignore() do
            if settings.bVerbose
                #println("calling update_inletQ_TotalQ")
            end
        end

        inletQ_TotalQ_local = nothing

        if nInletQ_BCs > 0
            inletQ_TotalQ_local = update_inletQ_TotalQ(inlet_discharges_current)
        end
    end

    Zygote.ignore() do
        if settings.bVerbose
            #@show typeof(inletQ_TotalQ_local)
            #@show inletQ_TotalQ_local
        end
    end

    # Process boundaries: update ghost cells values. Each boundary treatment function works on different part of Q_ghost. 
    Q_ghost_local = process_all_boundaries_2d(settings, Q, my_mesh_2D, boundary_conditions, ManningN_cells_local, zb_faces_local, swe_2D_constants, inletQ_Length, inletQ_TotalQ_local, exitH_WSE_local)

    Zygote.ignore() do
        if settings.bVerbose
            #@show typeof(Q_ghost_local)
            #@show size(Q_ghost_local)
            #@show Q_ghost_local
        end
    end

    # Loop through all cells to calculate the fluxes on faces
    updates = map(1:my_mesh_2D.numOfCells) do iCell           # .= is in-place mutation!
        cell_area = my_mesh_2D.cell_areas[iCell]

        # Initialize flux accumulation
        flux_sum = zeros(data_type, 3)

        for iFace in 1:my_mesh_2D.cellNodesCount[iCell]
            faceID = abs(my_mesh_2D.cellFacesList[iCell, :][iFace])
            left_cellID = iCell
            right_cellID = abs(my_mesh_2D.cellNeighbors_Dict[iCell][iFace])

            faceBoundaryID = my_mesh_2D.faceBoundaryID_Dict[faceID]
            face_normal = my_mesh_2D.cell_normals[iCell][iFace]

            if faceBoundaryID == 0  # internal face
                hL, huL, hvL = h[left_cellID], q_x[left_cellID], q_y[left_cellID]
                hR, huR, hvR = h[right_cellID], q_x[right_cellID], q_y[right_cellID]
            else  # boundary face
                hL, huL, hvL = h[left_cellID], q_x[left_cellID], q_y[left_cellID]
                hR = Q_ghost_local[right_cellID, 1]
                huR = Q_ghost_local[right_cellID, 2]
                hvR = Q_ghost_local[right_cellID, 3]
            end

            # Zygote.ignore() do
            #     if iCell == 1  # Print for first cell only
            #         println("Before Riemann solver:")
            #         @show hL, huL, hvL
            #         @show hR, huR, hvR
            #         @show face_normal
            #     end
            # end

            if RiemannSolver == "Roe"
                flux = Riemann_2D_Roe(settings, hL, huL, hvL, hR, huR, hvR, g, face_normal, hmin=h_small)
            elseif RiemannSolver == "HLL"
                error("HLL solver not implemented yet")
            elseif RiemannSolver == "HLLC"
                error("HLLC solver not implemented yet")
            else
                error("Wrong choice of RiemannSolver")
            end

            # Zygote.ignore() do
            #     if iCell == 1  # Print for first cell only
            #         println("After Riemann solver:")
            #         @show flux
            #     end
            # end

            # Accumulate flux contribution
            flux_sum = flux_sum .+ flux .* my_mesh_2D.face_lengths[faceID]
        end

        Zygote.ignore() do
            if iCell == -1
                println("flux_sum value = ", ForwardDiff.value.(flux_sum))
                println("flux_sum partials = ", ForwardDiff.partials.(flux_sum))
            end
        end

        # Source terms
        source_terms = if h[iCell] <= h_small
            [data_type(0.0),
                g * h[iCell] * S0[iCell, 1] * cell_area,
                g * h[iCell] * S0[iCell, 2] * cell_area]
        else
            u_temp = q_x[iCell] / max(h[iCell], h_small)
            v_temp = q_y[iCell] / max(h[iCell], h_small)
            u_mag = sqrt(u_temp^2 + v_temp^2 + sqrt(eps(data_type)))

            friction_x = g * ManningN_cells[iCell]^2 / (max(h[iCell], h_small))^(1.0 / 3.0) * u_mag * u_temp
            friction_y = g * ManningN_cells[iCell]^2 / (max(h[iCell], h_small))^(1.0 / 3.0) * u_mag * v_temp

            [zero(data_type),
                (g * h[iCell] * S0_local[iCell, 1] - friction_x) * cell_area,
                (g * h[iCell] * S0_local[iCell, 2] - friction_y) * cell_area]
        end

        if iCell == -1
            println("source_terms value = ", ForwardDiff.value.(source_terms))
            println("source_terms partials = ", ForwardDiff.partials.(source_terms))
        end

        if iCell == -1  # Print for first cell only
            Zygote.ignore() do
                if settings.bVerbose
                    @show flux_sum
                    @show source_terms
                    @show cell_area
                    @show (-flux_sum .+ source_terms) ./ cell_area
                end
            end
        end

        # Return the update for this cell (without in-place mutation)
        Array((-flux_sum .+ source_terms) ./ cell_area)

    end

    #convert updates to a 2D array: vcat organizes vectors as rows
    # Stacks vectors vertically, treating each vector as a row of the resulting matrix.
    # The transpose (') ensures the vectors are treated as rows when concatenated.
    #dQdt = vcat(updates'...)  #this breaks reverse mode AD in Zygote
    dQdt = [updates[i][j] for i in 1:length(updates), j in 1:length(updates[1])]
    # Ensure concrete type for output
    #dQdt = convert(typeof(Q), [updates[i][j] for i in 1:length(updates), j in 1:length(updates[1])])

    #dQdt = hcat(updates...)    #wrong dimensions
    #dQdt = reduce(hcat, updates) #wrong dimensions

    Zygote.ignore() do  
        if settings.bVerbose
            #@show size(updates)
            #@show updates
            #@show typeof(dQdt)
            #@show size(dQdt)
            #@show dQdt
        end
    end

    # Ensure output dQdt has same dimensions as Q
    @assert size(dQdt) == size(Q) "Dimension mismatch: dQdt $(size(dQdt)) ≠ Q $(size(Q))"

    return dQdt
end
