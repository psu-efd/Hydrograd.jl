
#define the Jacobian sparsity pattern
function define_sparsity_swe_2D(settings::ControlSettings, my_mesh_2D::mesh_2D)

    if settings.bVerbose
        println("define_sparsity_swe_2D")
    end

    jac_sparsity = nothing

    if (settings.bPerform_Forward_Simulation && settings.forward_settings.ode_solver_b_jac_sparsity) ||
       (settings.bPerform_Inversion && settings.inversion_settings.ode_solver_b_jac_sparsity) ||
       (settings.bPerform_Sensitivity_Analysis && settings.sensitivity_analysis_settings.ode_solver_b_jac_sparsity)

        # Populate the Jacobian sparsity pattern
        # Assume each variable depends on itself and its neighbors
        # we have 3 variables (h, q_x, q_y) for each cell
        jac_sparsity = spzeros(3 * my_mesh_2D.numOfCells, 3 * my_mesh_2D.numOfCells)

        for cellID in 1:my_mesh_2D.numOfCells
            # Self-dependence (diagonal entries)
            jac_sparsity[3*(cellID-1)+1, 3*(cellID-1)+1] = 1.0  # h -> h
            jac_sparsity[3*(cellID-1)+2, 3*(cellID-1)+2] = 1.0  # hu -> hu
            jac_sparsity[3*(cellID-1)+3, 3*(cellID-1)+3] = 1.0  # hv -> hv

            # Neighbor-dependence
            for neighbor in my_mesh_2D.cellNeighbors_Dict[cellID]

                #only need interior neighbors
                if neighbor > 0
                    jac_sparsity[3*(cellID-1)+1, 3*(neighbor-1)+1] = 1.0  # h -> h (neighbor)
                    jac_sparsity[3*(cellID-1)+2, 3*(neighbor-1)+2] = 1.0  # hu -> hu (neighbor)
                    jac_sparsity[3*(cellID-1)+3, 3*(neighbor-1)+3] = 1.0  # hv -> hv (neighbor)
                end
            end
        end
    end

    return jac_sparsity
end
