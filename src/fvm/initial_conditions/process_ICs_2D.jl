# process initial condition, such as setting up initial free surface, water depth, etc. 
# This file should be problem specific because each problem should have different ICs. 

# setup initial condition: wse, h, q_x, q_y
function setup_initial_condition!(settings, my_mesh_2D, nodeCoordinates, wse, zb_cell, h, q_x, q_y,swe_2D_constants, case_path, bPlot::Bool=false)

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
    else
        error("Invalid bPerform_Forward_Simulation, bPerform_Inversion, bPerform_Sensitivity_Analysis. No initial condition is to be setup.")
    end

    if settings.bVerbose
        println("setup_initial_condition!")
        println("initial_condition_options = ", initial_condition_options)
        println("initial_condition_constant_values = ", initial_condition_constant_values)
        println("initial_condition_values_from_file = ", initial_condition_values_from_file)
    end

    if initial_condition_options == "constant"
        #loop over cells
        for i in 1:my_mesh_2D.numOfCells
            h[i] = initial_condition_constant_values[1]
            q_x[i] = initial_condition_constant_values[2]
            q_y[i] = initial_condition_constant_values[3]
        end
    elseif initial_condition_options == "from_file"
        #make sure the length of the initial condition values is the same as the number of cells
        if length(initial_condition_values_from_file["h"]) != my_mesh_2D.numOfCells
            error("The length of the initial condition values of h (", length(initial_condition_values_from_file["h"]), ") is not the same as the number of cells (", my_mesh_2D.numOfCells, ").")
        end
        if length(initial_condition_values_from_file["q_x"]) != my_mesh_2D.numOfCells
            error("The length of the initial condition values of q_x (", length(initial_condition_values_from_file["q_x"]), ") is not the same as the number of cells (", my_mesh_2D.numOfCells, ").")
        end
        if length(initial_condition_values_from_file["q_y"]) != my_mesh_2D.numOfCells
            error("The length of the initial condition values of q_y (", length(initial_condition_values_from_file["q_y"]), ") is not the same as the number of cells (", my_mesh_2D.numOfCells, ").")
        end

        #copy the initial condition values
        copyto!(h, Float64.(initial_condition_values_from_file["h"]))
        copyto!(q_x, Float64.(initial_condition_values_from_file["q_x"]))
        copyto!(q_y, Float64.(initial_condition_values_from_file["q_y"]))
    end
   
    
    h[h.<0.0] .= swe_2D_constants.h_small  #ensure positivity of water depth h
    
    #update the free surface elevation again in case h has been clipped
    #wse = h + zb_cell
    copyto!(wse, h + zb_cell) 

    if settings.bVerbose
        println("initial conditions are setup.")
        println("wse = ", wse)
        println("h = ", h)
        println("q_x = ", q_x)
        println("q_y = ", q_y)
    end
  
    #optionally plot the ICs for checking 
    if bPlot
        vector_data = [] 
        vector_names = []
        
        scalar_data = [wse, h, q_x, q_y, zb_cell]
        scalar_names = ["wse", "h", "q_x", "q_y", "zb_cell"]
        
        file_path = joinpath(case_path, "initial_conditions.vtk" ) 
        export_to_vtk_2D(file_path, nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, 
                         scalar_data, scalar_names, vector_data, vector_names)    

        if settings.bVerbose
            println("initial conditions are saved to ", file_path)
        end
    end
    
end

#update ghost cells for wse, h, q_x, q_y
function setup_ghost_cells_initial_condition!(settings, my_mesh_2D, wse, h, q_x, q_y, wse_ghostCells, h_ghostCells, q_x_ghostCells, q_y_ghostCells)

    if settings.bVerbose
        println("setup_ghost_cells_initial_condition ...")
    end

    for iBoundaryFace in 1:my_mesh_2D.numOfAllBounaryFaces
        cellID_neighbor = my_mesh_2D.faceCells_Dict[my_mesh_2D.allBoundaryFacesIDs_List[iBoundaryFace]][1]
        
        wse_ghostCells[iBoundaryFace] = wse[cellID_neighbor]
        h_ghostCells[iBoundaryFace] = h[cellID_neighbor]
        q_x_ghostCells[iBoundaryFace] = q_x[cellID_neighbor]
        q_y_ghostCells[iBoundaryFace] = q_y[cellID_neighbor]
    end
    
    if settings.bVerbose
        println("Initital conditions for ghost cells are setup.")
        println("wse_ghostCells: ", wse_ghostCells)
        println("h_ghostCells: ", h_ghostCells)
        println("q_x_ghostCells: ", q_x_ghostCells)
        println("q_y_ghostCells: ", q_y_ghostCells)
    end
    
    
end 


