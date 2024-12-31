
#perform forward simulation for 2D SWE and save the results
function swe_2D_forward_simulation(settings, my_mesh_2D, swe_2D_constants, ode_f, Q0, params_array, 
            zb_cell_truth, ManningN_cell_truth, inlet_discharges_truth,
            total_water_volume, nodeCoordinates, zb_cells,
            save_path)

    #define the time for saving the results (for forward simulation, which may be different from the time for saving the results in inversion)
    dt_save = (swe_2D_constants.tspan[2] - swe_2D_constants.tspan[1]) / settings.forward_settings.nSave
    t_save = swe_2D_constants.tspan[1]:dt_save:swe_2D_constants.tspan[2]

    if settings.bVerbose
        println("swe_2D_forward_simulation")
        println("swe_2D_constants.tspan = ", swe_2D_constants.tspan)    
        println("swe_2D_constants.dt = ", swe_2D_constants.dt)
        println("dt_save = ", dt_save)
        println("t_save = ", t_save)
    end

    #define the ODE problem
    prob = ODEProblem(ode_f, Q0, swe_2D_constants.tspan, params_array)

    if settings.forward_settings.solver == "SciML"

        println("       with SciML solver ...")

        if settings.forward_settings.ode_solver == "Tsit5()"
            sol = solve(prob, Tsit5(), adaptive=settings.forward_settings.ode_solver_adaptive, dt=swe_2D_constants.dt, saveat=t_save)
        else
            sol = solve(prob, Tsit5(), adaptive=true, dt=swe_2D_constants.dt, saveat=t_save)
        end

        # #save the simulation solution results
        jldsave(joinpath(save_path, settings.forward_settings.save_file_name); sol)

        println(" Postprocessing the forward simulation results ...")
        
        # postprocess the simulation results
        postprocess_forward_simulation_results_swe_2D(settings, my_mesh_2D, total_water_volume, zb_cell_truth, ManningN_cell_truth, inlet_discharges_truth, zb_cells, nodeCoordinates, save_path)

    elseif settings.forward_settings.solver == "customized"    #My own ODE Solver
        #TODO: implement the customized ODE solver

        println("   Performing 2D SWE simulation with MyOwn solver ...")

        #sol = custom_ODE_solve(params_array, Q0, my_mesh_2D, swe_2D_constants.tspan, swe_2D_constants.dt)

        # #save the results
        # #save the simulation solution results
        #jldsave(joinpath(save_path, "simulation_solution.jld2"); sol)

        #swe_2D_save_results_custom(sol, total_water_volume, my_mesh_2D, zb_cells, save_path)
    else
        println("   Wrong solver choice. Supported solvers: SciML, customized. No forward simulation is performed.")
    end
end
