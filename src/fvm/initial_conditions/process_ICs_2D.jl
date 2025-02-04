# process initial condition, such as setting up initial free surface, water depth, etc. 
# This file should be problem specific because each problem should have different ICs. 

# setup initial condition: wse, wstill, h, hstill, xi, q_x, q_y
function setup_initial_condition(settings::ControlSettings, my_mesh_2D::mesh_2D, nodeCoordinates::Matrix{T}, wse::AbstractVector{T},
    wstill::AbstractVector{T}, h::AbstractVector{T}, hstill::AbstractVector{T}, xi::AbstractVector{T}, 
    q_x::AbstractVector{T}, q_y::AbstractVector{T},
    zb_cells::AbstractVector{T}, ManningN_cells::AbstractVector{T}, ks_cells::AbstractVector{T},
    swe_2D_constants::swe_2D_consts, case_path::String, bPlot::Bool=false) where {T<:Real}

    initial_condition_options = nothing
    initial_condition_constant_values = nothing
    initial_condition_values_from_file = nothing

    if settings.bPerform_Forward_Simulation
        initial_condition_options = settings.forward_settings.initial_condition_options
        initial_condition_constant_values = settings.forward_settings.initial_condition_constant_values
        initial_condition_values_from_file = settings.forward_settings.initial_condition_values_from_file
    elseif settings.bPerform_Inversion
        initial_condition_options = settings.inversion_settings.forward_simulation_initial_condition_options
        initial_condition_constant_values = settings.inversion_settings.forward_simulation_initial_condition_constant_values
        initial_condition_values_from_file = settings.inversion_settings.forward_simulation_initial_condition_values_from_file
    elseif settings.bPerform_Sensitivity_Analysis
        initial_condition_options = settings.sensitivity_analysis_settings.forward_simulation_initial_condition_options
        initial_condition_constant_values = settings.sensitivity_analysis_settings.forward_simulation_initial_condition_constant_values
        initial_condition_values_from_file = settings.sensitivity_analysis_settings.forward_simulation_initial_condition_values_from_file
    elseif settings.bPerform_UDE
        initial_condition_options = settings.UDE_settings.forward_simulation_initial_condition_options
        initial_condition_constant_values = settings.UDE_settings.forward_simulation_initial_condition_constant_values
        initial_condition_values_from_file = settings.UDE_settings.forward_simulation_initial_condition_values_from_file
    else
        error("Invalid bPerform_Forward_Simulation, bPerform_Inversion, bPerform_Sensitivity_Analysis, or bPerform_UDE. No initial condition is to be setup.")
    end

    #if settings.bVerbose
    println("setup_initial_condition")
    println("initial_condition_options = ", initial_condition_options)
    #println("initial_condition_constant_values = ", initial_condition_constant_values)
    #println("initial_condition_values_from_file = ", initial_condition_values_from_file)
    #if initial_condition_options == "from_file"
    #    println("initial_condition_values_from_file = ", length(initial_condition_values_from_file["wse"]))
    #    println("initial_condition_values_from_file = ", length(initial_condition_values_from_file["wstill"]))
    #    println("initial_condition_values_from_file = ", length(initial_condition_values_from_file["q_x"]))
    #    println("initial_condition_values_from_file = ", length(initial_condition_values_from_file["q_y"]))
    #end
    #end

    if initial_condition_options == "constant"
        wse_constant, wstill_constant, q_x_constant, q_y_constant = initial_condition_constant_values

        #loop over cells
        for i in 1:my_mesh_2D.numOfCells
            wstill[i] = wstill_constant
            wse[i] = wse_constant
            hstill[i] = wstill[i] - zb_cells[i]
            xi[i] = wse[i] - wstill[i]
            h[i] = wse[i] - zb_cells[i]

            q_x[i] = q_x_constant
            q_y[i] = q_y_constant
        end

    elseif initial_condition_options == "from_file"
        #make sure the length of the initial condition values is the same as the number of cells
        if length(initial_condition_values_from_file["wse"]) != my_mesh_2D.numOfCells
            error("The length of the initial condition values of wse (", length(initial_condition_values_from_file["wse"]), ") is not the same as the number of cells (", my_mesh_2D.numOfCells, ").")
        end
        if length(initial_condition_values_from_file["wstill"]) != my_mesh_2D.numOfCells
            error("The length of the initial condition values of wstill (", length(initial_condition_values_from_file["wstill"]), ") is not the same as the number of cells (", my_mesh_2D.numOfCells, ").")
        end
        if length(initial_condition_values_from_file["q_x"]) != my_mesh_2D.numOfCells
            error("The length of the initial condition values of q_x (", length(initial_condition_values_from_file["q_x"]), ") is not the same as the number of cells (", my_mesh_2D.numOfCells, ").")
        end
        if length(initial_condition_values_from_file["q_y"]) != my_mesh_2D.numOfCells
            error("The length of the initial condition values of q_y (", length(initial_condition_values_from_file["q_y"]), ") is not the same as the number of cells (", my_mesh_2D.numOfCells, ").")
        end

        #copy the initial condition values
        copyto!(wse, Float64.(initial_condition_values_from_file["wse"]))
        copyto!(wstill, Float64.(initial_condition_values_from_file["wstill"]))
        copyto!(q_x, Float64.(initial_condition_values_from_file["q_x"]))
        copyto!(q_y, Float64.(initial_condition_values_from_file["q_y"]))

        hstill .= wstill .- zb_cells
        xi .= wse .- wstill
        h .= wse .- zb_cells
    end   

    #make sure initial conditions make sense
    #if given wse is smaller than zb_cell, report an error
    if any(wse .< zb_cells)
        error("The initial condition for WSE is smaller than bed elevation (", wse, ").")
    end


    #if any of the h is smaller than h_small, report an error
    if any(h .< swe_2D_constants.h_small)
        error("The initial condition for water depth h from the given WSE and bed elevation is smaller than h_small (", swe_2D_constants.h_small, ").")
    end

    #h[h.<0.0] .= swe_2D_constants.h_small  #ensure positivity of water depth h

    if settings.bVerbose
        println("initial conditions are setup.")
        println("wse = ", wse)
        println("wstill = ", wstill)
        println("hstill = ", hstill)
        println("h = ", h)
        println("xi = ", xi)
        println("q_x = ", q_x)
        println("q_y = ", q_y)
    end

    #optionally plot the ICs for checking 
    if bPlot
        vector_data = []
        vector_names = []

        scalar_data = [wse, wstill, h, hstill, xi, q_x, q_y, zb_cells, ManningN_cells]
        scalar_names = ["wse", "wstill", "h", "hstill", "xi", "q_x", "q_y", "zb_cells", "ManningN_cells"]

        file_path = joinpath(case_path, "initial_conditions.vtk")
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount,
            "", "float", 0.0, scalar_data, scalar_names, vector_data, vector_names)

        if settings.bVerbose
            println("initial conditions are saved to ", file_path)
        end
    end

end

#update ghost cells for wse, wstill, h, hs, xi, q_x, q_y
function setup_ghost_cells_initial_condition(settings::ControlSettings, my_mesh_2D::mesh_2D, wse::AbstractVector{T}, 
    wstill::AbstractVector{T}, h::AbstractVector{T}, hstill::AbstractVector{T}, xi::AbstractVector{T}, 
    q_x::AbstractVector{T}, q_y::AbstractVector{T}, 
    wse_ghostCells::AbstractVector{T}, wstill_ghostCells::AbstractVector{T}, h_ghostCells::AbstractVector{T}, 
    hstill_ghostCells::AbstractVector{T}, xi_ghostCells::AbstractVector{T}, 
    q_x_ghostCells::AbstractVector{T}, q_y_ghostCells::AbstractVector{T}) where {T<:Real}

    if settings.bVerbose
        println("setup_ghost_cells_initial_condition ...")
    end

    for iBoundaryFace in 1:my_mesh_2D.numOfAllBounaryFaces
        cellID_neighbor = my_mesh_2D.faceCells_Dict[my_mesh_2D.allBoundaryFacesIDs_List[iBoundaryFace]][1]

        wse_ghostCells[iBoundaryFace] = wse[cellID_neighbor]
        wstill_ghostCells[iBoundaryFace] = wstill[cellID_neighbor]
        h_ghostCells[iBoundaryFace] = h[cellID_neighbor]
        hstill_ghostCells[iBoundaryFace] = hstill[cellID_neighbor]
        xi_ghostCells[iBoundaryFace] = xi[cellID_neighbor]
        q_x_ghostCells[iBoundaryFace] = q_x[cellID_neighbor]
        q_y_ghostCells[iBoundaryFace] = q_y[cellID_neighbor]
    end

    if settings.bVerbose
        println("Initital conditions for ghost cells are setup.")
        println("wse_ghostCells: ", wse_ghostCells)
        println("wstill_ghostCells: ", wstill_ghostCells)
        println("h_ghostCells: ", h_ghostCells)
        println("hstill_ghostCells: ", hstill_ghostCells)
        println("xi_ghostCells: ", xi_ghostCells)
        println("q_x_ghostCells: ", q_x_ghostCells)
        println("q_y_ghostCells: ", q_y_ghostCells)
    end


end


