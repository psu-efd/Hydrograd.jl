#semin-discretize the 2D SWE to convert it to an ODE system
#This file should be problem specific because we may invert different parameters 

# Problem-specific function to calculate the RHS of ODE: 
# dQdt (= dhdt dq_xdt dq_ydt)
# Q = [h, q_x, q_y] is the solution vector. q_x = h*u_x, q_y = h*u_y
# Q_ghost = [h_ghost, q_x_ghost, q_y_ghost] is the ghost cell values

using ForwardDiff
using Zygote


function swe_2d_rhs(Q, params_array, t, settings,
    my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants, ManningN_cells, ManningN_ghostCells, inletQ_Length, inletQ_TotalQ, exitH_WSE,
    zb_cells, zb_ghostCells, zb_faces, S0)

    data_type = eltype(Q)

    Zygote.ignore() do
        println("within swe_2D_rhs, t =", t)
        @show typeof(Q)
        @show typeof(params_array)
        @show typeof(params_array.ManningN_list_param)
        @show data_type
        @show params_array.ManningN_list_param
    end

    #return zeros(size(Q)) #zeros(data_type, size(Q))
    #return Q .* 1.0
    
    

    Zygote.ignore() do
        println("asserting data_type = ", data_type)
        @assert data_type <: Real "data_type must be a subtype of Real for AD compatibility"
    end


    # Mesh data
    #numOfCells = my_mesh_2D.numOfCells
    #maxNumOfCellFaces = my_mesh_2D.maxNumOfCellFaces

    g = swe_2D_constants.g
    h_small = swe_2D_constants.h_small
    RiemannSolver = swe_2D_constants.RiemannSolver

    #h = @view Q[:, 1]  # water depth
    #q_x = @view Q[:, 2]  # discharge in x
    #q_y = @view Q[:, 3]  # discharge in y
    h = Q[:, 1]  # water depth
    q_x = Q[:, 2]  # discharge in x
    q_y = Q[:, 3]  # discharge in y

    # Get parameter values
    # Handle both ComponentArray and regular Array
    if params_array isa ComponentArray
        zb_cells_current = Array(params_array.zb_cells_param)
        ManningN_list_current = params_array.ManningN_list_param
        inlet_discharges_current = Array(params_array.inlet_discharges_param)
    else
        # Assuming the order matches ComponentArray structure
        n_cells = my_mesh_2D.numOfCells
        n_materials = srh_all_Dict["srhmat_numOfMaterials"]
        zb_cells_current = params_array[1:n_cells]
        ManningN_list_current = params_array[n_cells+1:n_cells+n_materials]
        inlet_discharges_current = params_array[n_cells+n_materials+1:end]
    end

    Zygote.ignore() do
        @show typeof(params_array.zb_cells_param)
        @show typeof(params_array.ManningN_list_param)
        @show typeof(params_array.inlet_discharges_param)

        @show typeof(zb_cells_current)
        @show typeof(ManningN_list_current)
        @show typeof(inlet_discharges_current)

        @show zb_cells_current
        @show ManningN_list_current
        @show inlet_discharges_current
    end

    # Make the return value depend on parameters
    #return Q * 1.0 .+ reshape(zb_cells_current, :, 1) .* ones(1, 3)

    # For the case of inversion or sensitivity analysis, and if zb is an active parameter, 
    # we need to interpolate zb from cell to face and compute bed slope at cells
    if (settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis) &&
       "zb" in settings.inversion_settings.active_param_names

        zb_ghostCells, zb_faces, S0 = interploate_zb_from_cell_to_face_and_compute_S0(my_mesh_2D, zb_cells_current)
        #println("zb_ghostCells = ", zb_ghostCells.values)
        #println("S0 = ", S0)
    end

    Zygote.ignore() do
        @show typeof(zb_ghostCells)
        @show typeof(zb_faces)
        @show typeof(S0)
        @show size(S0)
        @show S0
        @show reshape(S0[:, 1], :, 1)
        @show reshape(S0[:, 1], :, 1) .* ones(1, 3)
    end

    # Make the return value depend on parameters
    #return Q * 1.0 .+ reshape(S0[:,1], :, 1) .* ones(1, 3)

    # try
    #     gradient_output = Zygote.gradient(new_ManningN_values -> begin
    #             ManningN_cells, ManningN_ghostCells = update_ManningN(my_mesh_2D, srh_all_Dict, new_ManningN_values)
    #             sum(ManningN_cells)
    #         end, ManningN_list_current)

    #     Zygote.ignore() do
    #         println("gradeint for update_ManningN successful")
    #         println("gradient_output for ManningN = ", gradient_output)
    #         @show gradient_output

    #         #throw(ErrorException("Stopping here for debugging"))
    #     end
    # catch e
    #     println("gradient for update_ManningN failed")
    #     @show e

    #     throw(ErrorException("Stopping here for debugging"))
    # end

    #For the case of inversion or sensitivity analysis, and if ManningN is an active parameter, 
    # we need to update ManningN at cells and ghost cells
    if (settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis) &&
       "ManningN" in settings.inversion_settings.active_param_names

       Zygote.ignore() do
            println("calling update_ManningN")
       end

        ManningN_cells, ManningN_ghostCells = update_ManningN(my_mesh_2D, srh_all_Dict, ManningN_list_current)
    end

    Zygote.ignore() do
        @show typeof(ManningN_cells)
        @show typeof(ManningN_ghostCells)
        @show size(ManningN_cells)
        @show ManningN_cells
        @show reshape(ManningN_cells, :, 1)
        @show reshape(ManningN_cells, :, 1) .* ones(1, 3)
    end

    # Make the return value depend on parameters
    #return Q * 1.0 .+ reshape(S0[:, 1], :, 1) .* ones(1, 3) .+ reshape(ManningN_cells, :, 1) .* ones(1, 3)
    #return Q * 1.0 .+ reshape(S0[:,1], :, 1) .* ones(1, 3) .+ (ManningN_list_current[1] + ManningN_list_current[2]) .* ones(size(Q))

    # try
    #     gradient_output = Zygote.gradient(new_inlet_discharge_values -> begin
    #                 inletQ_TotalQ = update_inletQ_TotalQ(new_inlet_discharge_values)
    #                 sum(inletQ_TotalQ)
    #             end, inlet_discharges_current)
    
    #     Zygote.ignore() do
    #             println("gradeint for update_inletQ_TotalQ successful")
    #             println("gradient_output for update_inletQ_TotalQ = ", gradient_output)
    #             @show gradient_output
    
    #             #throw(ErrorException("Stopping here for debugging"))
    #     end
    # catch e
    #         println("gradient for update_inletQ_TotalQ failed")
    #         @show e
    
    #         throw(ErrorException("Stopping here for debugging"))
    # end

    #For the case of inversion or sensitivity analysis, and if Q is an active parameter, 
    # we need to update inletQ_TotalQ based on the provided inlet_discharges_current
    if (settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis) &&
       "Q" in settings.inversion_settings.active_param_names

        inletQ_TotalQ = update_inletQ_TotalQ(inlet_discharges_current)
    end

    Zygote.ignore() do
        @show typeof(inletQ_TotalQ)
        @show inletQ_TotalQ
        @show sum(inletQ_TotalQ)
    end

    # Make the return value depend on parameters
    #return Q * 1.0 .+ reshape(S0[:, 1], :, 1) .* ones(1, 3) .+ reshape(ManningN_cells, :, 1) .* ones(1, 3)
    #return Q * 1.0 .+ reshape(S0[:, 1], :, 1) .* ones(1, 3) .+ reshape(ManningN_cells, :, 1) .* ones(1, 3) + inletQ_TotalQ[1] .* ones(size(Q))
    #return Q * 1.0 .+ reshape(S0[:,1], :, 1) .* ones(1, 3) .+ (ManningN_list_current[1] + ManningN_list_current[2]) .* ones(size(Q))

    # try
    #     gradient_output = Zygote.gradient((Q) -> sum(process_all_boundaries_2d(Q, my_mesh_2D, boundary_conditions, ManningN_cells, zb_faces, swe_2D_constants, inletQ_Length, inletQ_TotalQ, exitH_WSE)), Q)

    #     Zygote.ignore() do
    #         println("gradient_output for process_all_boundaries_2d successful")
    #         println("gradient_output for process_all_boundaries_2d = ", gradient_output)
    #         @show gradient_output
    #     end

    # catch e
    #     println("gradient for process_all_boundaries_2d failed")
    #     @show e

    #     throw(ErrorException("Stopping here for debugging"))
    # end

    # Process boundaries: update ghost cells values. Each boundary treatment function works on different part of Q_ghost. 
    Q_ghost = process_all_boundaries_2d(Q, my_mesh_2D, boundary_conditions, ManningN_cells, zb_faces, swe_2D_constants, inletQ_Length, inletQ_TotalQ, exitH_WSE)

    Zygote.ignore() do
        @show typeof(Q_ghost)
    end

    # Loop through all cells to calculate the fluxes on faces
    updates = collect(map(1:my_mesh_2D.numOfCells) do iCell           # .= is in-place mutation!
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
        source_terms = if h[iCell] <= h_small
            [data_type(0.0),
                g * h[iCell] * S0[iCell, 1] * cell_area,
                g * h[iCell] * S0[iCell, 2] * cell_area]
        else
            u_temp = q_x[iCell] / h[iCell]
            v_temp = q_y[iCell] / h[iCell]
            u_mag = max(sqrt(u_temp^2 + v_temp^2), sqrt(eps(data_type)))

            friction_x = g * ManningN_cells[iCell]^2 / h[iCell]^(1.0 / 3.0) * u_mag * u_temp
            friction_y = g * ManningN_cells[iCell]^2 / h[iCell]^(1.0 / 3.0) * u_mag * v_temp

            [zero(data_type),
                (g * h[iCell] * S0[iCell, 1] - friction_x) * cell_area,
                (g * h[iCell] * S0[iCell, 2] - friction_y) * cell_area]
        end

        if iCell == -1
            println("source_terms value = ", ForwardDiff.value.(source_terms))
            println("source_terms partials = ", ForwardDiff.partials.(source_terms))
        end

        # Return the update for this cell (without in-place mutation)
        Array((-flux_sum .+ source_terms) ./ cell_area)
    end)

    Zygote.ignore() do
        @show typeof(updates)
    end

    try
        #gradient_output_2 = Zygote.gradient(x -> sum(hcat(x...)), updates)   
        gradient_output_2 = Zygote.gradient(x -> sum(reduce(hcat, x)), updates)

        Zygote.ignore() do
            println("gradient_output for updates = ", gradient_output_2)
            @show gradient_output_2
        end

    catch e
        println("gradient(x -> sum(hcat(updates...)), updates) failed")
        @show e
        throw(error("stop here"))
    end



    #convert updates to a 2D array: vcat organizes vectors as rows
    # Stacks vectors vertically, treating each vector as a row of the resulting matrix.
    # The transpose (') ensures the vectors are treated as rows when concatenated.
    #dQdt = vcat(updates'...)  #this breaks reverse mode AD in Zygote
    #dQdt = [updates[i][j] for i in 1:length(updates), j in 1:length(updates[1])]
    #dQdt = hcat(updates...)
    dQdt = reduce(hcat, updates)

    Zygote.ignore() do
        @show typeof(dQdt)
    end

    #println("dQdt value = ", ForwardDiff.value.(dQdt))
    #println("dQdt partials = ", ForwardDiff.partials.(dQdt))

    #println("dQ/dpara = ", ForwardDiff.partials(Q))
    #throw("stop here")

    #return dQdt

    #return Q .* 1.0
    return Q * 1.0 .+ reshape(S0[:, 1], :, 1) .* ones(1, 3) .+ reshape(ManningN_cells, :, 1) .* ones(1, 3)
end
