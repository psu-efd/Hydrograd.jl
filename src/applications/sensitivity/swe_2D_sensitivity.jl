#perform sensitivity analysis for 2D SWE: Only the sensitivity of the final solution with respect to the active parameters is computed.
# For other kind of sensitivity analysis, a custom function should be implemented.

function swe_2D_sensitivity(ode_f, Q0, params_vector, swe_extra_params)

    #unpack the extra parameters
    active_param_name = swe_extra_params.active_param_name
    settings = swe_extra_params.settings
    my_mesh_2D = swe_extra_params.my_mesh_2D
    swe_2D_constants = swe_extra_params.swe_2D_constants
    nodeCoordinates = swe_extra_params.nodeCoordinates
    case_path = swe_extra_params.case_path
    zb_cells = swe_extra_params.zb_cells

    #get the still water surface elevation and still water depth
    wstill = swe_extra_params.wstill    
    hstill = swe_extra_params.hstill

    println("   Active parameter name = ", active_param_name)

    #define the time for saving the results (for sensitivity analysis, which may be different from the time for saving the results in forward simulation)
    dt_save = (swe_2D_constants.tspan[2] - swe_2D_constants.tspan[1]) / settings.sensitivity_analysis_settings.ode_solver_nSave
    t_save = swe_2D_constants.tspan[1]:dt_save:swe_2D_constants.tspan[2]

    if settings.bVerbose
        println("   swe_2D_sensitivity")
        println("   swe_2D_constants.tspan = ", swe_2D_constants.tspan)    
        println("   swe_2D_constants.dt = ", swe_2D_constants.dt)
        println("   dt_save = ", dt_save)
        println("   t_save = ", t_save)
    end

    function forward_simulation(θ)

        #define the ODE problem
        prob = ODEProblem(ode_f, Q0, swe_2D_constants.tspan, θ)

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
        pred_final = Array(pred)[:,end]

        Zygote.ignore() do
            #@show pred.retcode
            #@show typeof(pred)
            #@show size(pred)
            #@show pred

            #@show size(pred_final)

            #save the forward simulation results in a JSON file
            pred_array = Array(pred)
            #convert from ForwardDiff.Dual to Float64
            pred_array = map(x -> ForwardDiff.value(x), pred_array)
            
            #@show typeof(pred_array) 
            #@show size(pred_array)
            #@show pred_array            

            open(joinpath(case_path, "forward_simulation_results.json"), "w") do io
                JSON3.pretty(io, Dict("forward_simulation_results" => pred_array, "zb_cells" => zb_cells, "wstill" => wstill, "hstill" => hstill))
                println(io)
            end
            
        end

        return pred_final
    end

    #perform the sensitivity analysis (only ForwardDiff is supported for now)
    sensitivity = ForwardDiff.jacobian(forward_simulation, params_vector)

    #@show typeof(sensitivity)
    #@show size(sensitivity)
    
    #save the inversion results
    jldsave(joinpath(case_path, settings.sensitivity_analysis_settings.save_file_name); sensitivity)

    #also save the sensitivity results in a JSON file
    open(joinpath(case_path, "sensitivity_results.json"), "w") do io
        JSON3.pretty(io, Dict("sensitivity_results" => sensitivity, 
                              "parameter_name" => active_param_name,
                              "params_vector" => params_vector))
        println(io)
    end

    #process the inversion results
    println("   Post-processing sensitivity results ...")

    #process sensitivity results
    Hydrograd.postprocess_sensitivity_results_swe_2D(settings, swe_extra_params, my_mesh_2D, nodeCoordinates, params_vector, active_param_name, case_path)

end
