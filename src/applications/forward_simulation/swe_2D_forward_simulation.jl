
#perform forward simulation for 2D SWE and save the results
function swe_2D_forward_simulation(ode_f, Q0, params_vector, swe2d_extra_params, 
            zb_cell_truth, ManningN_zone_values_truth, inlet_discharges_truth
            )

    settings = swe2d_extra_params.settings
    swe_2D_constants = swe2d_extra_params.swe_2D_constants
    case_path = swe2d_extra_params.case_path

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
    prob = ODEProblem(ode_f, Q0, swe_2D_constants.tspan, params_vector)

    #debug AD correctness
    #debug start
    #params_vector_debug = [0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03]
    #debug_AD(ode_f, Q0, swe_2D_constants, params_vector_debug, swe2d_extra_params)
    #return
    #debug end

    if settings.forward_settings.solver == "SciML"

        println("       with SciML solver ...")

        if settings.forward_settings.ode_solver == "Tsit5()"
            sol = solve(prob, Tsit5(), adaptive=settings.forward_settings.ode_solver_adaptive, dt=swe_2D_constants.dt, saveat=t_save)
        else
            sol = solve(prob, Tsit5(), adaptive=true, dt=swe_2D_constants.dt, saveat=t_save)
        end

        # #save the simulation solution results
        save_ode_solution(sol, joinpath(case_path, settings.forward_settings.save_file_name))

        println("   Postprocessing the forward simulation results ...")
        
        # postprocess the simulation results
        postprocess_forward_simulation_results_swe_2D(swe2d_extra_params, zb_cell_truth, ManningN_zone_values_truth, inlet_discharges_truth)

    elseif settings.forward_settings.solver == "customized"    #My own ODE Solver
        #TODO: implement the customized ODE solver

        println("   Performing 2D SWE simulation with MyOwn solver ...")
        println("   This is not implemented yet.")

        #@show Q0
        
        #Profile.clear()
        #@profile sol = custom_ODE_solve(ode_f, Q0, params_vector, swe2d_extra_params)
        sol = custom_ODE_solve(ode_f, Q0, params_vector, swe2d_extra_params)
        #StatProfilerHTML.statprofilehtml()

        # #save the results
        # #save the simulation solution results
        jldsave(joinpath(swe2d_extra_params.case_path, "simulation_solution.jld2"); sol)

        swe_2D_save_results_custom(sol, swe2d_extra_params)
    else
        println("   Wrong solver choice. Supported solvers: SciML, customized. No forward simulation is performed.")
    end
end
