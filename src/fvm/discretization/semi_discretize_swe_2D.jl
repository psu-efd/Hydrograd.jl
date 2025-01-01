#semin-discretize the 2D SWE to convert it to an ODE system
#This file should be problem specific because we may invert different parameters 

# Problem-specific function to calculate the RHS of ODE: 
# dQdt (= dhdt dq_xdt dq_ydt)
# Q = [h, q_x, q_y] is the solution vector. q_x = h*u_x, q_y = h*u_y
# Q_ghost = [h_ghost, q_x_ghost, q_y_ghost] is the ghost cell values

# Arguments with "_passed" are passed to the function. They are only 
# used for forward simulations. For inversion and sensitivity analysis, 
# ManningN_cells, ManningN_ghostCells, inletQ_TotalQ, exitH_WSE, 
# zb_cells, zb_ghostCells, zb_faces, and S0 are all derived from params_vector.

#function swe_2d_rhs(Q, params_vector, active_range, param_ranges, t, settings,
#    my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants, inletQ_Length, exitH_WSE,
#    )

#function swe_2d_rhs(Q, params_vector, active_range, param_ranges, t, settings,
#    my_mesh_2D, srh_all_Dict, boundary_conditions, swe_2D_constants,
#    ManningN_cells_passed, ManningN_ghostCells_passed, inletQ_Length, inletQ_TotalQ_passed, exitH_WSE_passed,
#    zb_cells_passed, zb_ghostCells_passed, zb_faces_passed, S0_passed)

# @noinline function swe_2d_rhs(Q::Matrix{T1},
#     params_vector::Vector{T2},
#     active_range::UnitRange{Int},
#     param_ranges::Vector{UnitRange{Int}},
#     t::Real,
#     settings::ControlSettings,
#     my_mesh_2D::mesh_2D,
#     srh_all_Dict::Dict,
#     boundary_conditions::BoundaryConditions2D,
#     swe_2D_constants::swe_2D_consts,
#     ManningN_cells_passed::Vector{T3},
#     ManningN_ghostCells_passed::Vector{T4},
#     inletQ_Length_passed::Vector{Array{T5}},
#     inletQ_TotalQ_passed::Vector{T6},
#     exitH_WSE_passed::Vector{T7},
#     zb_cells_passed::Vector{T8},
#     zb_ghostCells_passed::Vector{T9},
#     zb_faces_passed::Vector{T10},
#     S0_passed::Matrix{T11}) where {T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}

using JuliaInterpreter

function swe_2d_rhs(Q::Matrix{T}, params_vector::Vector{T}, t::Float64, p_extra::Hydrograd.SWE2D_Extra_Parameters{T}) where {T}

    # Unpack the extra parameters
    active_param_name = p_extra.active_param_name
    settings = p_extra.settings
    my_mesh_2D = p_extra.my_mesh_2D
    srh_all_Dict = p_extra.srh_all_Dict
    boundary_conditions = p_extra.boundary_conditions
    swe_2D_constants = p_extra.swe_2D_constants
    ManningN_cells_passed = p_extra.ManningN_cells
    ManningN_ghostCells_passed = p_extra.ManningN_ghostCells
    inletQ_Length_passed = p_extra.inletQ_Length
    inletQ_TotalQ_passed = p_extra.inletQ_TotalQ
    exitH_WSE_passed = p_extra.exitH_WSE
    zb_cells_passed = p_extra.zb_cells
    zb_ghostCells_passed = p_extra.zb_ghostCells
    zb_faces_passed = p_extra.zb_faces
    S0_passed = p_extra.S0

    Zygote.ignore() do
        Base.@debug "swe_2d_rhs called" # This can help prevent inlining
    end

    # These are the arguments that are passed to the function. They are only 
    # used for forward simulations. For inversion and sensitivity analysis, 
    # ManningN_cells, ManningN_ghostCells, inletQ_TotalQ, exitH_WSE, 
    # zb_cells, zb_ghostCells, zb_faces, and S0 are all derived from params_vector (computed later in the function)
    ManningN_cells_local = deepcopy(ManningN_cells_passed)             #create a reference to ManningN_cells_passed. This makes Zygote happy (Don't know why)
    ManningN_ghostCells_local = deepcopy(ManningN_ghostCells_passed)
    inletQ_Length_local = deepcopy(inletQ_Length_passed)
    inletQ_TotalQ_local = deepcopy(inletQ_TotalQ_passed)
    exitH_WSE_local = deepcopy(exitH_WSE_passed)
    zb_cells_local = deepcopy(zb_cells_passed)                #not used for anything
    zb_ghostCells_local = deepcopy(zb_ghostCells_passed)
    zb_faces_local = deepcopy(zb_faces_passed)
    S0_local = deepcopy(S0_passed)

    g = swe_2D_constants.g
    h_small = swe_2D_constants.h_small
    RiemannSolver = swe_2D_constants.RiemannSolver

    #get the data type of Q
    data_type = eltype(Q)

    Zygote.ignore() do
        if settings.bVerbose
            #println("within swe_2D_rhs, t =", t)
            #println("asserting data_type = ", data_type)
        end

        @assert data_type <: Real "data_type must be a subtype of Real for AD compatibility"
    end

    Zygote.ignore() do
        if settings.bVerbose
            #@show typeof(params_vector)
            @show params_vector
        end

    end

    #return zeros(size(Q)) #zeros(data_type, size(Q))
    #return Q .* 1.0
   

    #h = @view Q[:, 1]  # water depth      #view of 2D array does not create a new array; only a reference
    #q_x = @view Q[:, 2]  # discharge in x
    #q_y = @view Q[:, 3]  # discharge in y
    h = copy(Q[:, 1])  # water depth        #slice of 2D array creates a new array (no need; waste of memory)
    q_x = copy(Q[:, 2])  # discharge in x
    q_y = copy(Q[:, 3])  # discharge in y

    # For the case of inversion or sensitivity analysis, and if zb is an active parameter, 
    # we need to interpolate zb from cell to face and compute bed slope at cells
    if ((settings.bPerform_Inversion ||settings.bPerform_Sensitivity_Analysis) && active_param_name == "zb")

        Zygote.ignore() do
            if settings.bVerbose
                println("calling interploate_zb_from_cell_to_face_and_compute_S0")
            end
        end

        zb_ghostCells_local, zb_faces_local, S0_local = interploate_zb_from_cell_to_face_and_compute_S0(my_mesh_2D, params_vector)
    end

    Zygote.ignore() do
        if settings.bVerbose
            @show typeof(zb_ghostCells_local)
            @show typeof(zb_faces_local)
            @show typeof(S0_local)
            @show size(S0_local)
            @show S0_local
        end
    end

    #For the case of inversion or sensitivity analysis, and if ManningN is an active parameter, 
    # we need to update ManningN at cells and ghost cells
    if (settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis) && active_param_name == "ManningN"

        Zygote.ignore() do
            if settings.bVerbose
                println("calling update_ManningN")
            end
        end

        ManningN_cells_local, ManningN_ghostCells_local = update_ManningN(my_mesh_2D, srh_all_Dict, params_vector)
    end

    Zygote.ignore() do
        if settings.bVerbose
            @show typeof(ManningN_cells_local)
            @show typeof(ManningN_ghostCells_local)
            @show size(ManningN_cells_local)
            @show ManningN_cells_local
        end
    end

    #For the case of inversion or sensitivity analysis, and if Q is an active parameter, 
    # we need to update inletQ_TotalQ based on the provided data
    if (settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis) && active_param_name == "Q"

        Zygote.ignore() do
            if settings.bVerbose
                println("updating inletQ_TotalQ")
            end
        end

        inletQ_TotalQ_local = params_vector        
    end

    Zygote.ignore() do
        if settings.bVerbose
            @show typeof(inletQ_TotalQ_local)
            @show inletQ_TotalQ_local
        end
    end

    Zygote.ignore() do
        if settings.bVerbose
            @show typeof(settings)
            @show typeof(Q)
            @show typeof(my_mesh_2D)
            @show typeof(boundary_conditions)
            @show typeof(ManningN_cells_local)
            @show typeof(zb_faces_local)
            @show typeof(swe_2D_constants)
            @show typeof(inletQ_Length_local)
            @show typeof(inletQ_TotalQ_local)
            @show typeof(exitH_WSE_local)
        end
    end

    #out = Zygote.@code_adjoint process_all_boundaries_2d(settings, h, q_x, q_y, my_mesh_2D, boundary_conditions, ManningN_cells_local, zb_faces_local, swe_2D_constants, 
    #                                          inletQ_Length_local, inletQ_TotalQ_local, exitH_WSE_local)
    #@code_warntype process_all_boundaries_2d(settings, h, q_x, q_y, my_mesh_2D, boundary_conditions, ManningN_cells_local, zb_faces_local, swe_2D_constants, 
    #                                          inletQ_Length_local, inletQ_TotalQ_local, exitH_WSE_local)
    #@show out

    # Process boundaries: update ghost cells values. Each boundary treatment function works on different part of Q_ghost. 
    h_ghost_local, q_x_ghost_local, q_y_ghost_local = process_all_boundaries_2d(settings, h, q_x, q_y, my_mesh_2D, boundary_conditions, 
                                              ManningN_cells_local, zb_cells_local, zb_faces_local, swe_2D_constants, 
                                              inletQ_Length_local, inletQ_TotalQ_local, exitH_WSE_local)

    Zygote.ignore() do
        if settings.bVerbose
            @show typeof(h_ghost_local)
            @show typeof(q_x_ghost_local)
            @show typeof(q_y_ghost_local)
            @show h_ghost_local
            @show q_x_ghost_local
            @show q_y_ghost_local
        end
    end

    # Loop through all cells to calculate the fluxes on faces
    updates = map(1:my_mesh_2D.numOfCells) do iCell           
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
                hR = h_ghost_local[right_cellID]
                huR = q_x_ghost_local[right_cellID]
                hvR = q_y_ghost_local[right_cellID]
            end

             Zygote.ignore() do
                 if iCell == -1  # Print for first cell only
                     println("Before Riemann solver:")
                     @show typeof(hL), typeof(huL), typeof(hvL)
                     @show hL, huL, hvL
                     @show typeof(hR), typeof(huR), typeof(hvR)
                     @show hR, huR, hvR
                    @show typeof(face_normal)
                    @show face_normal
                 end
             end

            if RiemannSolver == "Roe"
                flux = Riemann_2D_Roe(settings, hL, huL, hvL, hR, huR, hvR, g, face_normal, hmin=h_small)
            elseif RiemannSolver == "HLL"
                error("HLL solver not implemented yet")
            elseif RiemannSolver == "HLLC"
                error("HLLC solver not implemented yet")
            else
                error("Wrong choice of RiemannSolver")
            end

            Zygote.ignore() do
                if iCell == -1  # Print for first cell only
                    println("After Riemann solver:")
                    @show typeof(flux)
                    @show flux
                end
            end

            Zygote.ignore() do
                if iCell == -1  # Print for first cell only
                    println("before accumulating flux_sum")
                    @show typeof(my_mesh_2D.face_lengths)
                    @show my_mesh_2D.face_lengths
                    @show typeof(my_mesh_2D.face_lengths[faceID])
                    @show my_mesh_2D.face_lengths[faceID]
                    @show typeof(flux_sum)
                    @show flux_sum
                end
            end

            # Accumulate flux contribution
            flux_sum = flux_sum .+ flux .* my_mesh_2D.face_lengths[faceID]

            Zygote.ignore() do
                if settings.bVerbose
                    if iCell == -1  # Print for first cell only
                        println("after accumulating flux_sum")
                        @show typeof(flux_sum)
                        @show flux_sum
                    end
                end
            end
        end

        Zygote.ignore() do
            if settings.bVerbose
                if iCell == -1
                    @show typeof(flux_sum)
                    @show flux_sum
                end
            end
        end

        # Source terms
        source_terms = if h[iCell] <= h_small
            [zero(data_type),
                g * h[iCell] * S0_local[iCell, 1] * cell_area,
                g * h[iCell] * S0_local[iCell, 2] * cell_area]
        else
            u_temp = q_x[iCell] / max(h[iCell], h_small)
            v_temp = q_y[iCell] / max(h[iCell], h_small)
            u_mag = sqrt(u_temp^2 + v_temp^2 + sqrt(eps(data_type)))

            friction_x = g * ManningN_cells_local[iCell]^2 / (max(h[iCell], h_small))^(1.0 / 3.0) * u_mag * u_temp
            friction_y = g * ManningN_cells_local[iCell]^2 / (max(h[iCell], h_small))^(1.0 / 3.0) * u_mag * v_temp

            [zero(data_type),
                (g * h[iCell] * S0_local[iCell, 1] - friction_x) * cell_area,
                (g * h[iCell] * S0_local[iCell, 2] - friction_y) * cell_area]
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
            @show typeof(updates)
            @show size(updates)
            @show updates
            @show typeof(dQdt)
            @show size(dQdt)
            @show dQdt
        end
    end

    # Ensure output dQdt has same dimensions as Q
    @assert size(dQdt) == size(Q) "Dimension mismatch: dQdt $(size(dQdt)) ≠ Q $(size(Q))"

    return dQdt
end
