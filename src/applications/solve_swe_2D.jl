"""
    solve_swe_2D(control_file::String)

Solve 2D shallow water equations using settings from the specified control file.
The function supports forward simulation, inversion, sensitivity analysis, and UDE modes.

# Arguments
- `control_file::String`: Path to the JSON control file

# Returns
- `elapsed_seconds::Float64`: Total computation time in seconds
"""

function solve_swe_2D(control_file::String)

    #using Random
    #Random.seed!(1234)

    #hack for debugging
    #get the current directory
    #current_dir = pwd()
    #cd("case_path")
    
    # Add some warnings/info for gradient computation
    SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true

    #Timing 
    start_time = now()  # Current date and time

    #include the files to define the variables
    #get the part to the current file ("define_variables_swe_2D.jl" is in the same directory as the current file)
    current_file = @__FILE__
    current_dir = dirname(current_file)
    include(joinpath(current_dir, "define_variables_swe_2D.jl"))

    #print the banner
    Hydrograd.print_banner()

    #directory to read from/save to (the same directory as the current directory)
    #save_path = dirname(@__FILE__)
    case_path = pwd()

    # Read and parse control file
    println("Reading control file...")
    control_file = joinpath(case_path, "run_control.json")
    settings = Hydrograd.parse_control_file(control_file)

    #read data from SRH-2D hydro, geom, and material files; it aslo create the 2D mesh.
    srh_all_Dict = Hydrograd.process_SRH_2D_input(settings, case_path)

    #default unit system is SI
    unit_system = "SI"
    k_n = 1.0
    g = 9.81

    #define a swe_2D_constants object with some values from the control file
    swe_2D_constants = nothing
    if settings.time_settings.bUse_srhhydro_time_settings
        srhhydro_SimTime = srh_all_Dict["srhhydro_SimTime"]
        tStart = srhhydro_SimTime[1]
        tEnd = srhhydro_SimTime[2] * 3600 #convert hours to seconds
        dt = srhhydro_SimTime[3]

        grid_unit = srh_all_Dict["srhgeom_obj"].GridUnit        

        if grid_unit == "Meters"
            unit_system = "SI"
            k_n = 1.0
            g = 9.81
        elseif grid_unit == "Feet"
            unit_system = "EN"
            k_n = 1.486
            g = 32.2
        else
            error("Invalid grid unit: $grid_unit. Supported options: meter, ft.")
        end


        swe_2D_constants = Hydrograd.swe_2D_consts(
            unit_system,               #The choice of unit system, e.g., SI, EN
            k_n,               #unit conversion factor for Manning's equation
            g,               #gravity constant 
            tStart,              #current time 
            0.001,                #minimum time step size 
            dt,                   #time step size 
            0.4,                  #CFL number 
            tStart,               #start time for simulation
            tEnd,                 #end time for simulation 
            (tStart, tEnd),      #time span for simulation
            1.0e-3,               #dry bed water depth threshold (e.g., 1.0e-3)
            "Roe")                #The choice of Riemann solver, e.g., Roe, HLL, and HLLC                          

    else
        swe_2D_constants = Hydrograd.swe_2D_consts(
            unit_system,               #The choice of unit system, e.g., SI, EN
            k_n,               #unit conversion factor for Manning's equation
            g,               #gravity constant 
            settings.time_settings.tspan[1],              #current time 
            0.001,                #minimum time step size 
            settings.time_settings.dt,                   #time step size 
            0.4,                  #CFL number 
            settings.time_settings.tspan[1],  #start time for simulation
            settings.time_settings.tspan[2],  #end time for simulation 
            settings.time_settings.tspan,      #time span for simulation
            1.0e-3,               #dry bed water depth threshold (e.g., 1.0e-3)
            "Roe")                #The choice of Riemann solver, e.g., Roe, HLL, and HLLC                          
    end

    #Get the 2D mesh 
    my_mesh_2D = srh_all_Dict["my_mesh_2D"]

    # define solution variables at cells
    wse = zeros(my_mesh_2D.numOfCells)          #free surface elevation at cells 
    wstill = zeros(my_mesh_2D.numOfCells)      #still water elevation at cells 
    h = zeros(my_mesh_2D.numOfCells)            #water depth at cells 
    hstill = zeros(my_mesh_2D.numOfCells)      #still water depth at cells 
    xi = zeros(my_mesh_2D.numOfCells)          #free surface elevation above still water elevation at cells 
    q_x = zeros(my_mesh_2D.numOfCells)          #q_x=hu at cells 
    q_y = zeros(my_mesh_2D.numOfCells)          #q_y=hv at cells 

    # define solution variables at ghost cells
    wse_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #free surface elevation at ghost cells 
    wstill_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)      #still water elevation at ghost cells 
    hstill_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)      #still water depth at ghost cells 
    xi_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #free surface elevation above still water elevation at ghost cells 
    h_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)            #water depth at ghost cells 
    q_x_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #q_x=hu at ghost cells 
    q_y_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #q_y=hv at ghost cells 

    #mesh related variables which might be mutable (should not be in struct mesh_2D)
    #Get nodeCoordinates: Float64 2D array [numOfNodes, 3]
    nodeCoordinates = srh_all_Dict["srhgeom_obj"].nodeCoordinates

    #copy the bed elevation from nodeCoordinates to zb_nodes
    #zb_nodes = copy(nodeCoordinates[:, 3])

    #  setup bed elevation 
    #If performing inversion on zb, zero out the bed elevation in nodeCoordinates.
    if settings.bPerform_Inversion && settings.inversion_settings.active_param_names == ["zb"]
        nodeCoordinates[:, 3] .= 0.0
    end

    #setup bed elevation: compute zb at cell centers (zb_cells) from nodes, 
    #then interpolate zb from cell to face (zb_faces) and ghost cells (zb_ghostCells),
    #and compute the bed slope at cells (S0_cells) and faces (S0_faces). 
    #Note: If performing inversion on zb_cells, these values are not used in the inversion process. 
    #      Instead, they are updated in the inversion process.
    zb_cells, zb_ghostCells, zb_faces, S0_cells, S0_faces = Hydrograd.setup_bed(settings, my_mesh_2D, nodeCoordinates, case_path, true)

    #define the true bed elevation at cells and nodes
    zb_cells_truth = zeros(size(zb_cells))    

    #If performing forward simulation, make a copy of the bathymetry truth: zb_cell_truth; otherwise for 
    #inversion and sensitivity analysis, it is zero (its value will be loaded from the forward simulation result in the inversion).
    if settings.bPerform_Forward_Simulation
        zb_cells_truth = deepcopy(zb_cells)
    end

    #get the true Manning's n and inlet discharges
    srhhydro_ManningsN_Dict = srh_all_Dict["srhhydro_ManningsN"]

    #define the true Manning's n values for all material zones
    ManningN_zone_values_truth = zeros(length(srhhydro_ManningsN_Dict))

    #If performing forward simulation, make a copy of the Manning's n truth: ManningN__zone_values_truth; otherwise for 
    #inversion and sensitivity analysis, it is zero (its value will be loaded from the forward simulation result in the inversion).
    if settings.bPerform_Forward_Simulation
        ManningN_zone_values_truth = Float64[srhhydro_ManningsN_Dict[i] for i in 0:(length(srhhydro_ManningsN_Dict)-1)]  #0 is default material ID.
    end

    #get the true inlet discharges (could be nothing if no inletQ_BCs)
    srhhydro_inletQ_Dict = srh_all_Dict["srhhydro_IQParams"]

    #define the true inlet discharges
    inlet_discharges_truth = Vector{Float64}()

    if length(srhhydro_inletQ_Dict) > 0
        #@show srhhydro_inletQ_Dict
        inlet_discharges_truth = [parse(Float64, srhhydro_inletQ_Dict[k][1]) for k in keys(srhhydro_inletQ_Dict)]
    end

    #print the true values. The true parameter values are computed regardless of whether performing forward simulation or inversion
    if settings.bVerbose
        println("True zb_cells = ", zb_cells_truth)
        println("True zb_nodes = ", zb_nodes_truth)
        println("True Manning's n values = ", ManningN_zone_values_truth)
        println("True inlet discharges = ", inlet_discharges_truth)
    end

    #preprocess: create a 1D array for model parameter for 2D shallow water equations
    #params_vector: the 1D array of the active parameter (zb_cells_truth, ManningN_zone_values_truth, or inlet_discharges_truth)
    #active_param_name: the name of the active parameter (zb, ManningN, or Q)
    params_vector, active_param_name = Hydrograd.setup_model_parameters_2D(settings, my_mesh_2D, srh_all_Dict, zb_cells_truth, ManningN_zone_values_truth, inlet_discharges_truth)

    #Initial setup of Manning's n at cells and ghost cells using the SRH-2D data 
    #If performing forward simulation and the ManningN_option is "variable_as_function_of_h", ManningN_cells will be updated later in the forward simulation process.
    #If performing inversion on Manning's n, ManningN_cells will be updated later in the inversion process.
    #If performing UDE and UDE_choice is ManningN_h, ManningN_cells will be updated later in the UDE process.
    ManningN_cells = Hydrograd.setup_ManningN(settings, my_mesh_2D, srh_all_Dict)

    #setup initial condition for wse, wstill, h, hstill, xi, q_x, q_y at cells
    Hydrograd.setup_initial_condition(settings, my_mesh_2D, nodeCoordinates, wse, wstill, h, hstill, xi, q_x, q_y, zb_cells, ManningN_cells,
                                      swe_2D_constants, case_path, true)

    #setup initial condition for wse, h, q_x, q_y at ghost cells
    Hydrograd.setup_ghost_cells_initial_condition(settings, my_mesh_2D, wse, wstill, h, hstill, xi, q_x, q_y, 
                                      wse_ghostCells, wstill_ghostCells, h_ghostCells, hstill_ghostCells, xi_ghostCells, q_x_ghostCells, q_y_ghostCells)

    #update the dry/wet flag at cells and the flags for whether a cell is adjacent to dry land and high dry land
    b_dry_wet, b_Adjacent_to_dry_land, b_Adjacent_to_high_dry_land = Hydrograd.process_dry_wet_flags(my_mesh_2D, h, zb_cells, swe_2D_constants) 

    #@show typeof(b_dry_wet)
    #@show typeof(b_Adjacent_to_dry_land)
    #@show typeof(b_Adjacent_to_high_dry_land)

    #create and preprocess boundary conditions: boundary_conditions only contains the static information of the boundaries.
    boundary_conditions, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_Length, inletQ_TotalA, inletQ_DryWet,
    exitH_WSE, exitH_H, exitH_A, wall_H, wall_A, symm_H, symm_A = Hydrograd.initialize_boundary_conditions_2D(settings, srh_all_Dict, nodeCoordinates)

    #set up initial condition for for solution state variables for ODE solver
    Q0 = vcat(xi, q_x, q_y)   #Q0 is a 1D array

    # Define the UDE model (if not performing UDE, ude_model, ude_model_params, ude_model_state will not be used)
    ude_model, ude_model_params, ude_model_state = Hydrograd.create_NN_model(settings)

    #define whether to use in-place ODE solver
    bInPlaceODE = false #default is false. Use in-place ODE solver if the sensealg supports it.
    if (settings.bPerform_Inversion && (settings.inversion_settings.ode_solver_sensealg == "AutoEnzyme()" ||
                                        settings.inversion_settings.ode_solver_sensealg == "ForwardSensitivity()" ||
                                        settings.inversion_settings.ode_solver_sensealg == "ForwardDiffSensitivity()")) #|| settings.bPerform_Forward_Simulation
        bInPlaceODE = true
    end

    #@show bInPlaceODE

    # Create the extra parameters struct
    swe_extra_params = SWE2D_Extra_Parameters(
        case_path,
        active_param_name,
        settings,
        bInPlaceODE,
        my_mesh_2D,
        nodeCoordinates,
        srh_all_Dict,
        boundary_conditions,
        swe_2D_constants,
        ManningN_cells,
        inletQ_Length,
        inletQ_TotalQ,
        exitH_WSE,
        wse,
        wstill,
        h,
        hstill,
        wse_ghostCells,
        wstill_ghostCells,
        h_ghostCells,
        hstill_ghostCells,
        zb_cells,
        zb_ghostCells,
        zb_faces,
        S0_cells,
        S0_faces,
        b_dry_wet,
        b_Adjacent_to_dry_land,
        b_Adjacent_to_high_dry_land,
        ude_model,
        ude_model_params,
        ude_model_state
    )

    # Create the ODEFunction with the extra parameters struct passed in.
    #This is the out-of-place version of the ODE function.
    ode_f_out_of_place = ODEFunction((u, p, t) -> begin

            if settings.bVerbose
                #@show typeof(u)
                #@show typeof(p)
                #@show typeof(t)
                #@show typeof(swe_extra_params)
            end

            #dummpy du (not used)
            du = zeros(size(u))

            swe_2d_rhs(du, u, p, t, swe_extra_params)
        end; jac_prototype=jac_sparsity)

    #This is the in-place version of the ODE function.
    ode_f_in_place = ODEFunction((du, u, p, t) -> begin

            if settings.bVerbose
                #@show typeof(u)
                #@show typeof(p)
                #@show typeof(t)
                #@show typeof(swe_extra_params)
            end

            swe_2d_rhs(du, u, p, t, swe_extra_params)
        end; jac_prototype=jac_sparsity)


    #########################
    #Forward Simulation part#
    #########################

    if settings.bPerform_Forward_Simulation
        println("Forward simulation (2D SWE) ...")

        #perform forward simulation (out-of-place version)
        #Profile.clear()
        #@profile
        # For forward simulation, we always use out-of-place ODE solver (it seems to be slightly faster).
        if swe_extra_params.bInPlaceODE
            println("   Using in-place ODE solver for forward simulation ...")
            Hydrograd.swe_2D_forward_simulation(ode_f_in_place, Q0, params_vector, swe_extra_params,
                zb_cells_truth, ManningN_zone_values_truth, inlet_discharges_truth)
        else
            println("   Using out-of-place ODE solver for forward simulation ...")
            Hydrograd.swe_2D_forward_simulation(ode_f_out_of_place, Q0, params_vector, swe_extra_params,
                zb_cells_truth, ManningN_zone_values_truth, inlet_discharges_truth)
        end

        #StatProfilerHTML.statprofilehtml()

        #Profile.print()

        #ProfileView.view()

        #ProfileView.save("profile.svg")  # Save as SVG

        #readline()
    end


    #########################
    #Inversion part         # 
    #########################

    if settings.bPerform_Inversion

        println("Parameter inversion ...")

        #perform inversion
        #@code_warntype 
        if swe_extra_params.bInPlaceODE
            Hydrograd.swe_2D_inversion(ode_f_in_place, Q0, params_vector, swe_extra_params)
        else
            Hydrograd.swe_2D_inversion(ode_f_out_of_place, Q0, params_vector, swe_extra_params)
        end

    end


    #########################   
    #Sensitivity part       # 
    #########################

    if settings.bPerform_Sensitivity_Analysis

        println("Sensitivity analysis ...")

        #perform sensitivity analysis (out-of-place version)
        Hydrograd.swe_2D_sensitivity(ode_f_out_of_place, Q0, params_vector, swe_extra_params)

    end


    #########################
    #UDE part               # 
    #########################

    if settings.bPerform_UDE

        if settings.UDE_settings.UDE_mode == "training"
            println("UDE training ...")
        elseif settings.UDE_settings.UDE_mode == "inference"
            println("UDE inference ...")
        else
            error("Invalid UDE mode: $(settings.UDE_settings.UDE_mode). Supported options: training, inference.")
        end

        #perform UDE
        Hydrograd.swe_2D_UDE(ode_f, Q0, params_vector, swe_extra_params)
    end

    #Timing 
    end_time = now()  # Current date and time
    elapsed_time = end_time - start_time
    elapsed_seconds = Millisecond(elapsed_time).value / 1000
    println("Elapsed time in seconds: $elapsed_seconds")

    #restore the current directory
    #cd(current_dir)

    println("All done!")

    return elapsed_seconds

end
