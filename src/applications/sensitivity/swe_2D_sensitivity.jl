#perform sensitivity analysis for 2D SWE: Only the sensitivity of the final solution with respect to the active parameters is computed.

function swe_2D_sensitivity(settings, my_mesh_2D, swe_2D_constants, ode_f, Q0, params_array, active_range, param_ranges,
    nodeCoordinates, zb_cells, case_path)

    #define the time for saving the results (for sensitivity analysis, which may be different from the time for saving the results in forward simulation)
    dt_save = (swe_2D_constants.tspan[2] - swe_2D_constants.tspan[1]) / settings.sensitivity_analysis_settings.ode_solver_nSave
    t_save = swe_2D_constants.tspan[1]:dt_save:swe_2D_constants.tspan[2]

    if settings.bVerbose
        println("swe_2D_sensitivity")
        println("swe_2D_constants.tspan = ", swe_2D_constants.tspan)    
        println("swe_2D_constants.dt = ", swe_2D_constants.dt)
        println("dt_save = ", dt_save)
        println("t_save = ", t_save)
    end

    function forward_simulation(θ)

        #create the full parameter array
        p_new = [i ∈ active_range ? θ[i-active_range[1]+1] : convert(eltype(θ), params_array[i]) for i in 1:length(params_array)]

        prob = ODEProblem(ode_f, Q0, swe_2D_constants.tspan, p_new)

        # Solve the ODE (forward pass)
        if settings.sensitivity_analysis_settings.ode_solver == "Tsit5()"
            pred = solve(prob, Tsit5(), 
                adaptive=settings.sensitivity_analysis_settings.ode_solver_adaptive, 
                dt=swe_2D_constants.dt, saveat=t_save, 
                sensealg=SensitivityADPassThrough()) #working
        else
            #pred = solve(prob, Tsit5(), adaptive=true, dt=dt, saveat=t_save)
            error("Not implemented yet")
        end

        #get the final solution
        sol = pred[:,:,end]

        Zygote.ignore() do
            #@show pred.retcode
            #@show typeof(pred)
            #@show size(pred)
            #@show pred
        end

        return sol
    end


    #get the active sensitivity parameters: params_array
    θ = params_array[active_range]  #get a copy of the subarray of params_array corresponding to the active parameters

    #perform the sensitivity analysis
    sol = ForwardDiff.jacobian(forward_simulation, θ)
    
    #save the inversion results
    jldsave(joinpath(case_path, settings.sensitivity_analysis_settings.save_file_name); sol)

    #process the inversion results
    println("   Post-processing sensitivity results ...")

    #process sensitivity results
    Hydrograd.postprocess_sensitivity_results_swe_2D(settings, my_mesh_2D, nodeCoordinates, param_ranges, case_path)

end
