#semin-discretize the 2D SWE to convert it to an ODE system
#This file should be problem specific because we may invert different parameters 

# Problem-specific function to calculate the RHS of ODE: 
# dQdt (= dhdt dq_xdt dq_ydt)
# Q = [h, q_x, q_y] is the solution vector. q_x = h*u_x, q_y = h*u_y
# Q_ghost = [h_ghost, q_x_ghost, q_y_ghost] is the ghost cell values

using ForwardDiff
using Zygote


function swe_2d_rhs(Q, params_array, t, bPerform_Forward_Simulation, bPerform_Inversion, bPerform_Sensitivity_Analysis, 
    my_mesh_2D, boundary_conditions, swe_2D_constants, ManningN_cells, ManningN_ghostCells, inletQ_Length, inletQ_TotalQ, exitH_WSE,
    zb_cells, zb_ghostCells, zb_faces, S0, 
    active_params_names)

    #Zygote.ignore() do  
    #    println("within swe_2D_rhs, t =", t)
    #end

    # Mesh data
    #numOfCells = my_mesh_2D.numOfCells
    #maxNumOfCellFaces = my_mesh_2D.maxNumOfCellFaces

    g = swe_2D_constants.g
    h_small = swe_2D_constants.h_small
    RiemannSolver = swe_2D_constants.RiemannSolver

    h = @view Q[:, 1]  # water depth
    q_x = @view Q[:, 2]  # discharge in x
    q_y = @view Q[:, 3]  # discharge in y

    # Get parameter values
    zb_cells_current = params_array.zb_cells_param  
    ManningN_list_current = params_array.ManningN_list_param
    inlet_discharges_current = params_array.inlet_discharges_param

    # For the case of inversion or sensitivity analysis, and if zb is an active parameter, 
    # we need to interpolate zb from cell to face and compute bed slope at cells
    if (bPerform_Inversion || bPerform_Sensitivity_Analysis) && "zb" in active_params_names
        zb_ghostCells, zb_faces, S0 = interploate_zb_from_cell_to_face_and_compute_S0(my_mesh_2D, zb_cells_current)
        #println("zb_ghostCells = ", zb_ghostCells.values)
        #println("S0 = ", S0)
    end

    #For the case of inversion or sensitivity analysis, and if ManningN is an active parameter, 
    # we need to update ManningN at cells and ghost cells
    if (bPerform_Inversion || bPerform_Sensitivity_Analysis) && "ManningN" in active_params_names
        ManningN_cells, ManningN_ghostCells = update_ManningN(my_mesh_2D, ManningN_list_current)
    end

    #For the case of inversion or sensitivity analysis, and if Q is an active parameter, 
    # we need to update inletQ_TotalQ based on the provided inlet_discharges_current
    if (bPerform_Inversion || bPerform_Sensitivity_Analysis) && "Q" in active_params_names
        inletQ_TotalQ = update_inletQ_TotalQ(inlet_discharges_current)
    end

    # Process boundaries: update ghost cells values. Each boundary treatment function works on different part of Q_ghost. 
    Q_ghost = process_all_boundaries_2d(Q, my_mesh_2D, boundary_conditions, ManningN_cells, zb_faces, swe_2D_constants, inletQ_Length, inletQ_TotalQ, exitH_WSE)

    # Loop through all cells to calculate the fluxes on faces
    updates = map(1:my_mesh_2D.numOfCells) do iCell           # .= is in-place mutation!
        cell_area = my_mesh_2D.cell_areas[iCell]

        # Initialize flux accumulation
        flux_sum = zeros(eltype(Q), 3)

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
            flux_sum = flux_sum .+ flux .* my_mesh_2D.face_lengths[faceID]
        end

        if iCell == -1
            println("flux_sum value = ", ForwardDiff.value.(flux_sum))
            println("flux_sum partials = ", ForwardDiff.partials.(flux_sum))
        end

        # Source terms
        # source_terms = zeros(eltype(Q), 3)
        # if h[iCell] <= h_small
        #     source_terms .= [eltype(Q)(0.0),
        #         g * h[iCell] * S0[iCell, 1] * cell_area,
        #         g * h[iCell] * S0[iCell, 2] * cell_area]
        # else
        #     u_temp = q_x[iCell] / h[iCell]
        #     v_temp = q_y[iCell] / h[iCell]
        #     #u_mag = sqrt(u_temp^2 + v_temp^2)
        #     u_mag = max(sqrt(u_temp^2 + v_temp^2), sqrt(eps(eltype(u_temp))))

        #     friction_x = g * ManningN_cells[iCell]^2 / h[iCell]^(1.0 / 3.0) * u_mag * u_temp
        #     friction_y = g * ManningN_cells[iCell]^2 / h[iCell]^(1.0 / 3.0) * u_mag * v_temp

        #     source_terms .= [zero(eltype(Q)),
        #         (g * h[iCell] * S0[iCell, 1] - friction_x) * cell_area,
        #         (g * h[iCell] * S0[iCell, 2] - friction_y) * cell_area]
        # end

        source_terms = if h[iCell] <= h_small
            [eltype(Q)(0.0),
             g * h[iCell] * S0[iCell, 1] * cell_area,
             g * h[iCell] * S0[iCell, 2] * cell_area]
        else
            u_temp = q_x[iCell] / h[iCell]
            v_temp = q_y[iCell] / h[iCell]
            u_mag = max(sqrt(u_temp^2 + v_temp^2), sqrt(eps(eltype(u_temp))))
        
            friction_x = g * ManningN_cells[iCell]^2 / h[iCell]^(1.0 / 3.0) * u_mag * u_temp
            friction_y = g * ManningN_cells[iCell]^2 / h[iCell]^(1.0 / 3.0) * u_mag * v_temp
        
            [zero(eltype(Q)),
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
    #dQdt = vcat(updates'...)  #this breaks reverse mode AD in Zygote
    dQdt = [updates[i][j] for i in 1:length(updates), j in 1:length(updates[1])]

    #println("dQdt value = ", ForwardDiff.value.(dQdt))
    #println("dQdt partials = ", ForwardDiff.partials.(dQdt))

    #println("dQ/dpara = ", ForwardDiff.partials(Q))
    #throw("stop here")

    return dQdt
end
