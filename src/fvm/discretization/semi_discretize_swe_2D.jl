#semin-discretize the 2D SWE to convert it to an ODE system

# Problem-specific function to calculate the RHS of ODE: 
# dQdt (= dxidt dq_xdt dq_ydt)
# Q = [xi, q_x, q_y] is the solution vector. q_x = h*u_x, q_y = h*u_y

# Notes on the arguments and params_vector:
# For forward simulations: params_vector is not used at all. Model parameters (ManningN, zb, inletQ, etc.) are passed within extra parameters (p_extra).
# For inversion and sensitivity analysis: params_vector is the trainable parameters (ManningN, zb, or inletQ.) that are optimized or whose sensitivities are computed. In 
#          addition to the trainable parameters, other parameters are passed within extra parameters (p_extra).
# For UDE: params_vector is the trainable NN parameters. For ManningN_h option, the output of UDE(params_vector) gives the ManningN values for a given h. The ManningN 
#          values within p_extra are used.  
#          For flow_resistance option, the output of UDE(params_vector) directly gives the flow resistance values for a given (h, u_x, u_y; no need of ManningN at all).
#          ude_model, ude_model_params, and ude_model_state are passed within p_extra.

#using JuliaInterpreter

function swe_2d_rhs(dQdt::AbstractVector{T1}, Q::AbstractVector{T2}, params_vector::AbstractVector{T3}, 
    t::Float64, p_extra::Hydrograd.SWE2D_Extra_Parameters{T4})::AbstractVector{promote_type(T1, T2, T3, T4)} where {T1,T2,T3,T4}

    Zygote.ignore() do
        #@show t
    end

    # Unpack the extra parameters
    active_param_name = p_extra.active_param_name
    settings = p_extra.settings
    my_mesh_2D = p_extra.my_mesh_2D
    srh_all_Dict = p_extra.srh_all_Dict
    boundary_conditions = p_extra.boundary_conditions
    swe_2D_constants = p_extra.swe_2D_constants

    #unpack the dry/wet flags
    b_dry_wet = p_extra.b_dry_wet
    b_Adjacent_to_dry_land = p_extra.b_Adjacent_to_dry_land
    b_Adjacent_to_high_dry_land = p_extra.b_Adjacent_to_high_dry_land

    # These are the arguments that are passed to the function. They may be only 
    # used for forward simulations. For inversion and sensitivity analysis, 
    # ManningN_cells, inletQ_TotalQ, exitH_WSE, or
    # zb_cells, zb_ghostCells, zb_faces, and S0 are all derived from params_vector (computed later in the function)
    ManningN_cells_local = p_extra.ManningN_cells          #get a reference to these vectors
    inletQ_Length_local = p_extra.inletQ_Length
    inletQ_TotalQ_local = p_extra.inletQ_TotalQ
    exitH_WSE_local = p_extra.exitH_WSE
    zb_cells_local = p_extra.zb_cells
    zb_ghostCells_local = p_extra.zb_ghostCells
    zb_faces_local = p_extra.zb_faces
    S0_cells_local = p_extra.S0_cells
    S0_faces_local = p_extra.S0_faces

    hstill_local = p_extra.hstill
    wstill_local = p_extra.wstill
    
    hstill_ghost_local = p_extra.hstill_ghostCells
    wstill_ghost_local = p_extra.wstill_ghostCells

    #other variables from swe_2D_constants
    g = swe_2D_constants.g
    k_n = swe_2D_constants.k_n
    h_small = swe_2D_constants.h_small
    RiemannSolver = swe_2D_constants.RiemannSolver

    # Use promote_type to ensure consistent types for computations
    data_type = promote_type(T1, T2, T3, T4)

    Zygote.ignore() do
        if settings.bVerbose
            #println("within swe_2D_rhs, t =", t)
            #println("asserting data_type = ", data_type)
        end

        @assert data_type <: Real "data_type must be a subtype of Real for AD compatibility"
    end

    Zygote.ignore() do
        if settings.bVerbose
            if eltype(params_vector) <: ForwardDiff.Dual
                println("params_vector contains dual numbers")
                @show ForwardDiff.value.(params_vector)  # Get values from all dual numbers
            else
                println("params_vector contains regular numbers")
                @show params_vector
            end
        end
    end

    xi = @view Q[1:my_mesh_2D.numOfCells]  # water depth      #view of 2D array does not create a new array; only a reference
    q_x = @view Q[my_mesh_2D.numOfCells+1:2*my_mesh_2D.numOfCells]  # discharge in x
    q_y = @view Q[2*my_mesh_2D.numOfCells+1:3*my_mesh_2D.numOfCells]  # discharge in y

    #@show xi 
    #@show q_x
    #@show q_y

    h = xi .+ hstill_local

    # Replace in-place mutations with functional style
    h = ifelse.(h .<= h_small, h_small, h)
    q_x = ifelse.(h .<= h_small, zero(data_type), q_x)
    q_y = ifelse.(h .<= h_small, zero(data_type), q_y)

    #@show h
    #@show q_x
    #@show q_y

    # For the case of inversion or sensitivity analysis, and if zb is the active parameter, 
    # we need to update bed data from zb_cells
    if ((settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis) && active_param_name == "zb")

        Zygote.ignore() do
            if settings.bVerbose
                println("calling update_bed_data")
            end
        end

        #zb_cells_local points to params_vector, which is the zb parameter vector to be inverted
        zb_cells_local = params_vector

        zb_ghostCells_local, zb_faces_local, S0_cells_local, S0_faces_local = update_bed_data(my_mesh_2D, params_vector)        
    end

    Zygote.ignore() do
        #if settings.bVerbose
            #@show typeof(zb_ghostCells_local)
            #@show typeof(zb_faces_local)
            #@show typeof(S0_cells_local)
            #@show size(S0_cells_local)
            #@show S0_cells_local[1]
        #end
    end

    #For the case of forward simulation: if ManningN_option is constant, ManningN_cells_local is already updated in the preprocess step (no need to update here).
    #If ManningN_option is variable_as_function_of_h, ManningN_cells_local is updated here.
    if settings.bPerform_Forward_Simulation && settings.forward_settings.ManningN_option == "variable_as_function_of_h"
        #ManningN_cells_local is a new allocated array, so it is not the same as p_extra.ManningN_cells as the update. 
        ManningN_cells_local = update_ManningN_forward_simulation(h, settings)

        #For the case of inversion or sensitivity analysis, and if ManningN is the active parameter, 
        # we need to update ManningN at cells and ghost cells
    elseif (settings.bPerform_Inversion || settings.bPerform_Sensitivity_Analysis) && active_param_name == "ManningN"

        Zygote.ignore() do
            if settings.bVerbose
                println("calling update_ManningN_inversion_sensitivity_analysis")
            end
        end

        ManningN_cells_local = update_ManningN_inversion_sensitivity_analysis(my_mesh_2D, srh_all_Dict, params_vector)

        #For the case of UDE, we need to update ManningN or flow resistance based on the UDE model
        #In this case, params_vector is the trainable NN parameters
    elseif settings.bPerform_UDE && settings.UDE_settings.UDE_choice == "ManningN_h"
        #ManningN_cells_local = update_ManningN_UDE(h, p_extra.ude_model, p_extra.ude_model_params, p_extra.ude_model_state, my_mesh_2D)

        #@show typeof(Float64.(settings.UDE_settings.UDE_NN_config["h_bounds"]))
        #@show Float64.(settings.UDE_settings.UDE_NN_config["h_bounds"])

        ManningN_cells_local = update_ManningN_UDE(h, p_extra.ude_model, params_vector, p_extra.ude_model_state, Float64.(settings.UDE_settings.UDE_NN_config["h_bounds"]), my_mesh_2D.numOfCells)
    end

    Zygote.ignore() do
        if settings.bVerbose
            #@show typeof(ManningN_cells_local)
            #@show size(ManningN_cells_local)
            #@show ManningN_cells_local
        end
    end

    #For the case of inversion or sensitivity analysis, and if Q is the active parameter, 
    # we need to update inletQ_TotalQ based on the provided params_vector
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
            #@show typeof(inletQ_TotalQ_local)
            #@show inletQ_TotalQ_local
            #@show typeof(zb_cells_local)
            #@show typeof(zb_faces_local)
            
        end
    end

    #update the dry/wet flags
    b_dry_wet, b_Adjacent_to_dry_land, b_Adjacent_to_high_dry_land = process_dry_wet_flags(my_mesh_2D, h, zb_cells_local, swe_2D_constants)

    # Process boundaries: update ghost cells values. Each boundary treatment function works on different part of h_ghost, q_x_ghost, and q_y_ghost. 
    h_ghost_local, q_x_ghost_local, q_y_ghost_local = process_all_boundaries_2d(settings, h, q_x, q_y, my_mesh_2D, boundary_conditions,
        ManningN_cells_local, zb_cells_local, zb_faces_local, swe_2D_constants,
        inletQ_Length_local, inletQ_TotalQ_local, exitH_WSE_local)

    # Update ghost cells for xi
    xi_ghost_local = h_ghost_local .- hstill_ghost_local

    Zygote.ignore() do
        if settings.bVerbose
            #@show typeof(h_ghost_local)
            #@show typeof(q_x_ghost_local)
            #@show typeof(q_y_ghost_local)
            #@show h_ghost_local
            #@show q_x_ghost_local
            #@show q_y_ghost_local
        end
    end

    #compute the contribution of inviscid fluxes
    updates_inviscid = compute_inviscid_fluxes(settings, xi, hstill_local, h, q_x, q_y, xi_ghost_local, 
                            hstill_ghost_local, h_ghost_local, q_x_ghost_local, q_y_ghost_local, 
                            ManningN_cells_local, 
                            zb_cells_local, zb_ghostCells_local, zb_faces_local, S0_faces_local, my_mesh_2D, g, RiemannSolver, h_small, data_type)

    #@show typeof(updates_inviscid)
    #@show size(updates_inviscid)
    #@show updates_inviscid
    
    #@show typeof(updates_inviscid)
    #@show size(updates_inviscid)
    #@show updates_inviscid
    
    #compute the contribution of source terms
    updates_source = compute_source_terms(settings, my_mesh_2D, xi, h, q_x, q_y, S0_cells_local, ManningN_cells_local, 
                              params_vector, p_extra, g, k_n, h_small)

    #combine inviscid and source terms
    if p_extra.bInPlaceODE    
        # In-place: use broadcast assignment to update existing array
        #@. dQdt = updates_inviscid + updates_source  # This creates an intermediate array and then assigns it to dQdt.
        dQdt .= updates_inviscid .+ updates_source  # This is more efficient because it avoids the intermediate array creation. (not the case; even slower)
    #    return dQdt  # Explicit return for clarity
    else    
        # Out-of-place: create new array
        #@show size(updates_inviscid)
        #@show size(updates_source)
        #@show (updates_inviscid .+ updates_source)[7], (updates_inviscid .+ updates_source)[my_mesh_2D.numOfCells+7], 
        #       (updates_inviscid .+ updates_source)[2*my_mesh_2D.numOfCells+7]

        return updates_inviscid .+ updates_source
    end

    #@show typeof(Q)
    #@show size(Q)
    #@show typeof(dQdt)
    #@show size(dQdt)
    #@show dQdt

    # Ensure output dQdt has same dimensions as Q
    #@assert size(dQdt) == size(Q) "Dimension mismatch: dQdt $(size(dQdt)) ≠ Q $(size(Q))"

    #return dQdt
end

#function to compute the inviscid fluxes
#serial version
function compute_inviscid_fluxes(settings, xi, hstill, h, q_x, q_y, xi_ghost, hstill_ghost, h_ghost, q_x_ghost, q_y_ghost, ManningN_cells, 
    zb_cells_local, zb_ghostCells_local, zb_faces_local, S0_faces_local, 
    my_mesh_2D, g, RiemannSolver, h_small, data_type)

    #compute the inviscid fluxes for each cell (in the order of cells)
    updates_inviscid_cells = [
        let
            cell_area = my_mesh_2D.cell_areas[iCell]

            # Initialize flux accumulation
            flux_sum = zeros(data_type, 3)

            for iFace in 1:my_mesh_2D.cellNodesCount[iCell]
                faceID = my_mesh_2D.cellFacesList[iCell, :][iFace]
                left_cellID = iCell
                right_cellID = my_mesh_2D.cellNeighbors_Dict[iCell][iFace]

                faceBoundaryID = my_mesh_2D.faceBoundaryID_Dict[faceID]
                face_normal = my_mesh_2D.cell_normals[iCell][iFace]

                if !my_mesh_2D.bFace_is_boundary[faceID]  # internal face
                    #hL, huL, hvL = h[left_cellID]+zb_cells_local[left_cellID]-zb_faces_local[faceID], q_x[left_cellID], q_y[left_cellID]
                    #hR, huR, hvR = h[right_cellID]+zb_cells_local[right_cellID]-zb_faces_local[faceID], q_x[right_cellID], q_y[right_cellID]

                    xiL, hstillL, hL, huL, hvL = xi[left_cellID], hstill[left_cellID], h[left_cellID], q_x[left_cellID], q_y[left_cellID]
                    xiR, hstillR, hR, huR, hvR = xi[right_cellID], hstill[right_cellID], h[right_cellID], q_x[right_cellID], q_y[right_cellID]

                    zb_L = zb_cells_local[left_cellID]
                    zb_R = zb_cells_local[right_cellID]
                else  # boundary face
                    #hL, huL, hvL = h[left_cellID]+zb_cells_local[left_cellID]-zb_faces_local[faceID], q_x[left_cellID], q_y[left_cellID]
                    #hR = h_ghost[right_cellID]+zb_ghostCells_local[right_cellID]-zb_faces_local[faceID]

                    xiL = xi[left_cellID]
                    hstillL = hstill[left_cellID]
                    hL = h[left_cellID]
                    huL = q_x[left_cellID]
                    hvL = q_y[left_cellID]

                    xiR = xi_ghost[right_cellID]
                    hstillR = hstill_ghost[right_cellID]
                    hR = h_ghost[right_cellID]
                    huR = q_x_ghost[right_cellID]
                    hvR = q_y_ghost[right_cellID]

                    zb_L = zb_cells_local[left_cellID]
                    zb_R = zb_ghostCells_local[right_cellID]
                end

                #get face bed elevation
                zb_face = zb_faces_local[faceID]

                Zygote.ignore() do
                    if iCell == -8  
                        println("Before Riemann solver:")
                        @show iCell, iFace, left_cellID, right_cellID, my_mesh_2D.bFace_is_boundary[faceID]
                        #@show typeof(hL), typeof(huL), typeof(hvL)
                        @show hL, huL, hvL
                        #@show typeof(hR), typeof(huR), typeof(hvR)
                        @show hR, huR, hvR
                        #@show typeof(face_normal)
                        @show face_normal
                        
                    end
                end

                #left cell is the current cell
                S0_on_face = S0_faces_local[left_cellID][iFace]

                #Compute the inviscid fluxes using the Riemann solver. In fact, 
                #the return flux also contains the bed slope term contribution
                if RiemannSolver == "Roe"
                    flux = Riemann_2D_Roe(iCell, settings, xiL, hstillL, hL, huL, hvL, zb_L, 
                                                           xiR, hstillR, hR, huR, hvR, zb_R, 
                                   zb_face, S0_on_face, g, face_normal, hmin=h_small)
                elseif RiemannSolver == "HLL"
                    error("HLL solver not implemented yet")
                elseif RiemannSolver == "HLLC"
                    error("HLLC solver not implemented yet")
                else
                    error("Wrong choice of RiemannSolver")
                end

                Zygote.ignore() do
                    if iCell == -8  
                        println("After Riemann solver:")
                        @show iCell, iFace, left_cellID, right_cellID
                        #@show typeof(flux)
                        #@show size(flux)
                        @show flux
                    end
                end

                Zygote.ignore() do
                    if iCell == -7  # Print for first cell only
                        println("before accumulating flux_sum")
                        #@show typeof(my_mesh_2D.face_lengths)
                        #@show my_mesh_2D.face_lengths
                        #@show typeof(my_mesh_2D.face_lengths[faceID])
                        @show my_mesh_2D.face_lengths[faceID]
                        #@show typeof(flux_sum)
                        @show flux_sum
                    end
                end

                # Accumulate flux contribution
                flux_sum = flux_sum .+ flux .* my_mesh_2D.face_lengths[faceID]

                Zygote.ignore() do
                    #if settings.bVerbose
                        if iCell == -7  # Print for first cell only
                            println("after accumulating flux_sum")
                            #@show iCell, iFace, left_cellID, right_cellID
                            #@show typeof(flux_sum)
                            #@show size(flux_sum)
                            @show flux_sum
                        end
                    #end
                end
            end

            Zygote.ignore() do
                #if settings.bVerbose
                    if iCell == -8
                        println("after accumulating flux_sum")
                        @show iCell
                        #@show typeof(flux_sum)
                        #@show size(flux_sum)
                        @show -flux_sum ./ cell_area
                    end
                #end
            end

            # Return the update for this cell (without in-place mutation)            
            -flux_sum ./ cell_area

        end
        for iCell in 1:my_mesh_2D.numOfCells
    ]

    #convert updates_inviscid_cells to a 1D array in the order (h1, h2, h3, q_x1, q_x2, q_x3, q_y1, q_y2, q_y3)
    updates_inviscid = rearrange_vector_of_vectors(updates_inviscid_cells)

    Zygote.ignore() do
        #if settings.bVerbose
            #@show typeof(updates_inviscid)
            #@show size(updates_inviscid)
            #@show updates_inviscid
        #end
    end

    return updates_inviscid

end

#function to rearrange a vector of vectors into a 1D array
function rearrange_vector_of_vectors(data::Vector{Vector{T}}) where T
    n = length(data)  # Number of elements in the outer vector
    m = length(data[1])  # Number of components in each tuple

    # Separate components into individual arrays
    components = [getindex.(data, i) for i in 1:m]

    # Concatenate the components into a single 1D array
    return vcat(components...)
end

#function to compute the source terms
function compute_source_terms(settings::ControlSettings, my_mesh_2D::mesh_2D, xi::AbstractVector{T1}, h::AbstractVector{T1}, q_x::AbstractVector{T1}, q_y::AbstractVector{T1},
    S0::Matrix{T2}, ManningN_cells::Vector{T3}, params_vector::AbstractVector{T4}, p_extra::SWE2D_Extra_Parameters{T5}, g::Float64, k_n::Float64, h_small::Float64) where {T1,T2,T3,T4,T5}

    data_type = promote_type(T1, T2, T3, T4, T5)

    #unpack the dry/wet flags
    b_Adjacent_to_high_dry_land = p_extra.b_Adjacent_to_high_dry_land

    #compute the friction (flow resistance) terms
    friction_x, friction_y = compute_friction_terms(settings, h, q_x, q_y, ManningN_cells, params_vector, p_extra, my_mesh_2D, g, k_n, h_small)

    #compute the momentum source terms

    #create a boolean array to check if h is less than or equal to h_small
    above_small_h = h .> h_small

    # remove the intermediate arrays
    #source_x = below_small_h .* g .* h .* S0[:, 1] .+ .~below_small_h .* (g .* h .* S0[:, 1] .- friction_x)
    #source_y = below_small_h .* g .* h .* S0[:, 2] .+ .~below_small_h .* (g .* h .* S0[:, 2] .- friction_y)

    #source terms considering the dry/wet flags: if a wet cell is adjacent to high dry land, the bed slope part of source terms is zero.
    #source_x = above_small_h .* (g .* h .* S0[:, 1] .* .~b_Adjacent_to_high_dry_land .- friction_x)
    #source_y = above_small_h .* (g .* h .* S0[:, 2] .* .~b_Adjacent_to_high_dry_land .- friction_y)

    #source terms considering the bed slope term contribution
    source_x = above_small_h .* (g .* xi .* S0[:, 1] .- friction_x)
    source_y = above_small_h .* (g .* xi .* S0[:, 2] .- friction_y)
    
    # combine them into a single 1D array
    updates_source = vcat(zeros(data_type, length(h)), source_x, source_y)

    Zygote.ignore() do
        if settings.bVerbose
            
                println("source terms, iCell = 7")
                @show above_small_h[7]
                @show h[7], q_x[7], q_y[7]
                @show S0[7, 1], S0[7, 2]
                @show ManningN_cells[7]
                #@show params_vector[7]
                @show friction_x[7], friction_y[7]
                @show source_x[7], source_y[7]
                @show updates_source[7], updates_source[my_mesh_2D.numOfCells+7], updates_source[2*my_mesh_2D.numOfCells+7]
            
        end
    end

    Zygote.ignore() do
        #if settings.bVerbose
            #@show typeof(updates_source)
            #@show size(updates_source)
            #@show updates_source
        #end
    end

    return updates_source
end

#function to compute the friction (flow resistance) terms
function compute_friction_terms(settings::ControlSettings, h::AbstractVector{T1}, q_x::AbstractVector{T1}, q_y::AbstractVector{T1}, 
    ManningN_cells::Vector{T2}, params_vector::Union{AbstractVector{T3}, Nothing}, p_extra::SWE2D_Extra_Parameters{T4}, 
    my_mesh_2D::mesh_2D, g::Float64, k_n::Float64, h_small::Float64) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real} 

    #initialize friction terms
    #friction_x = zeros(eltype(h), my_mesh_2D.numOfCells)
    #friction_y = zeros(eltype(h), my_mesh_2D.numOfCells)

    #if performing UDE and UDE_choice is FlowResistance, compute the friction terms from the UDE model
    if settings.bPerform_UDE && settings.UDE_settings.UDE_choice == "FlowResistance"
        friction_magnitudes = update_FlowResistance_UDE(
            h, q_x, q_y,
            p_extra.ude_model,
            params_vector,
            p_extra.ude_model_state,
            Float64.(settings.UDE_settings.UDE_NN_config["h_bounds"]),
            Float64.(settings.UDE_settings.UDE_NN_config["velocity_magnitude_bounds"]),
            my_mesh_2D.numOfCells
        )

        mag_q = smooth_sqrt.(q_x.^2 .+ q_y.^2)

        #compute the friction terms
        friction_x = [friction_magnitudes[i] * q_x[i] / mag_q[i] for i in eachindex(h)]
        friction_y = [friction_magnitudes[i] * q_y[i] / mag_q[i] for i in eachindex(h)]

        #original code: check small h (no longer needed; friction terms are zero for dry cells anyway)
        #if the cell is dry or the magnitude of the velocity is too small, set the friction terms to zero
        #friction_x = [h[i] < h_small || mag_q[i] < 1e-6 ? zero(eltype(h)) :
        #              friction_magnitudes[i] * q_x[i] / mag_q[i]
        #              for i in eachindex(h)]
        #friction_y = [h[i] < h_small || mag_q[i] < 1e-6 ? zero(eltype(h)) :
        #              friction_magnitudes[i] * q_y[i] / mag_q[i]
        #              for i in eachindex(h)]
    else  # Just use the Manning's formula
        # Compute friction terms 
        friction_x = [g * ManningN_cells[i]^2/ k_n^2 / (h[i]+h_small)^(7.0 / 3.0) * smooth_sqrt(q_x[i]^2 + q_y[i]^2) * q_x[i]
                      for i in eachindex(h)]
        friction_y = [g * ManningN_cells[i]^2/ k_n^2 / (h[i]+h_small)^(7.0 / 3.0) * smooth_sqrt(q_x[i]^2 + q_y[i]^2) * q_y[i]
                      for i in eachindex(h)]
    
        #original code: check small h (no longer needed; friction terms are zero for dry cells anyway)
        #friction_x = [h[i] < h_small ? zero(T1) :
        #              g * ManningN_cells[i]^2/ k_n^2 / h[i]^(7.0 / 3.0) * smooth_sqrt(q_x[i]^2 + q_y[i]^2) * q_x[i]
        #              for i in eachindex(h)]
        #friction_y = [h[i] < h_small ? zero(T1) :
        #              g * ManningN_cells[i]^2 / k_n^2 / h[i]^(7.0 / 3.0) * smooth_sqrt(q_x[i]^2 + q_y[i]^2) * q_y[i]
        #              for i in eachindex(h)]
    end

    return friction_x, friction_y
end

#combine functions for AD debugging
# function combined_functions(settings, h, q_x, q_y, h_ghost, q_x_ghost, q_y_ghost, my_mesh_2D, S0, ManningN_cells, g, RiemannSolver, h_small, data_type)

#     updates_inviscid = compute_inviscid_fluxes(settings, h, q_x, q_y, h_ghost, q_x_ghost, q_y_ghost, ManningN_cells, my_mesh_2D, g, RiemannSolver, h_small, data_type)
#     updates_source = compute_source_terms(settings, my_mesh_2D, h, q_x, q_y, S0, ManningN_cells, g, h_small, data_type)

#     updates = updates_inviscid .+ updates_source

#     return updates

# end 
