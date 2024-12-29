# Custom ODE solver for SWE_2D


# Finite Volume Update
function custom_ODE_update_cells(Q, para, t, swe_2D_constants, my_mesh_2D)

    #dQdt = zeros(eltype(para), my_mesh_2D.numOfCells, 3)
    
    # compute dQdt 
    #swe_2d_ode!(dQdt, Q, para, t)
    #dQdt = swe_2d_ode(Q, para, t)
    dQdt = 0.0
    
     # Vectorized update of Q
    Q_new = Q + swe_2D_constants.dt * dQdt
    
    # Create mask for dry cells
    dry_mask = Q_new[:, 1] .< swe_2D_constants.h_small
    
    # Create new array with all updates at once
    Q_new = hcat(
        ifelse.(dry_mask, swe_2D_constants.h_small, Q_new[:, 1]),
        ifelse.(dry_mask, zero(eltype(para)), Q_new[:, 2]),
        ifelse.(dry_mask, zero(eltype(para)), Q_new[:, 3])
    )
   
    # Update WSE using broadcasting
    #eta = Q_new[:, 1] + zb_cells
   
    return Q_new
                
end

# Main Solver Function (equivalent to solve(prob, Euler(), adaptive=false, p=para, dt=dt, saveat=t_save))
function custom_ODE_solve(para, Q0, my_mesh_2D, swe_2D_constants)
    t_start = swe_2D_constants.tspan[1]
    t_end = swe_2D_constants.tspan[2]
    time_steps = t_start:swe_2D_constants.dt:t_end

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
        Q = custom_ODE_update_cells(Q, para, t, swe_2D_constants, my_mesh_2D)

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
