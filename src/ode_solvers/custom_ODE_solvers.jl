# Custom ODE solver for SWE_2D


# Finite Volume Update
function custom_ODE_update_cells(ode_f, Q, params_vector, t, dt, swe2d_extra_params)

    #dQdt = zeros(eltype(para), my_mesh_2D.numOfCells, 3)
    
    # compute dQdt 
    #swe_2d_ode!(dQdt, Q, para, t)
    dQdt = ode_f(Q, params_vector, t)
    
     # Vectorized update of Q
    Q_new = Q + dt * dQdt
    
    # Create mask for dry cells
    dry_mask = Q_new[:, 1] .< swe2d_extra_params.swe_2D_constants.h_small
    
    # Create new array with all updates at once
    Q_new = hcat(
        ifelse.(dry_mask, swe2d_extra_params.swe_2D_constants.h_small, Q_new[:, 1]),
        ifelse.(dry_mask, zero(eltype(params_vector)), Q_new[:, 2]),
        ifelse.(dry_mask, zero(eltype(params_vector)), Q_new[:, 3])
    )
   
    # Update WSE using broadcasting
    #eta = Q_new[:, 1] + zb_cells
   
    return Q_new
                
end

# Main Solver Function (equivalent to solve(prob, Euler(), adaptive=false, p=para, dt=dt, saveat=t_save))
function custom_ODE_solve(ode_f, Q0, params_vector, swe2d_extra_params)
    t_start = swe2d_extra_params.swe_2D_constants.tspan[1]
    t_end = swe2d_extra_params.swe_2D_constants.tspan[2]
    dt = swe2d_extra_params.swe_2D_constants.dt
    time_steps = t_start:dt:t_end

    nTimeSteps = length(time_steps)

    #save frequency (every how many time steps)
    save_freq = 1

    #how many saves 
    nSaves = nTimeSteps ÷ save_freq

    #define the array to store the solution
    #sol = zeros(eltype(para), my_mesh_2D.numOfCells, 3, nSaves)

    #list of saved solutions
    solution_list = []

    #initialize the solution
    Q = Q0

    iStep = 1
    
    for t in time_steps
        Zygote.ignore() do
            println("time = ", t)

            if iStep == 5000
                print("Here.")
            end
        end
        
        #update Q
        Q = custom_ODE_update_cells(ode_f, Q, params_vector, t, dt, swe2d_extra_params)

        #save the solution
        if (iStep+1) % save_freq == 0
            solution_list = [solution_list..., copy(Q)]  # Create new list instead of mutating
        end
                        
        iStep += 1
    end

    # Convert the list of solutions to a 3D array using cat
    sol = cat(solution_list..., dims=3)

    # sol is a 3D array with dimensions [cells, variables, timesteps]
    # Convert to tuple of 2D arrays (one per timestep) to make like a solution from DifferentialEquations.jl
    #sol_tuples = tuple(eachslice(sol, dims=3)...)

    return sol
    
end
