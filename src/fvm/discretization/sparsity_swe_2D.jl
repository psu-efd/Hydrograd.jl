
#define the Jacobian sparsity pattern
function define_sparsity_swe_2D(settings::ControlSettings, my_mesh_2D::mesh_2D)

    if settings.bVerbose
        println("define_sparsity_swe_2D")
    end

    jac_sparsity = nothing

    if (settings.bPerform_Forward_Simulation && settings.forward_settings.ode_solver_b_jac_sparsity) ||
       (settings.bPerform_Inversion && settings.inversion_settings.ode_solver_b_jac_sparsity) ||
       (settings.bPerform_Sensitivity_Analysis && settings.sensitivity_analysis_settings.ode_solver_b_jac_sparsity)

        numOfCells = my_mesh_2D.numOfCells

        # Populate the Jacobian sparsity pattern
        # Assume each variable depends on itself and its neighbors
        # we have 3 variables (h, q_x, q_y) for each cell
        jac_sparsity = spzeros(3 * numOfCells, 3 * numOfCells)


        for cellID in 1:numOfCells
            # Self-dependence (diagonal entries)
            jac_sparsity[cellID, cellID] = 1.0  # h -> h
            jac_sparsity[cellID+numOfCells, cellID+numOfCells] = 1.0  # hu -> hu
            jac_sparsity[cellID+2*numOfCells, cellID+2*numOfCells] = 1.0  # hv -> hv

            cell_faces = my_mesh_2D.cellFacesList[cellID, :]
            cellNeighbors = my_mesh_2D.cellNeighbors_Dict[cellID]

            # Neighbor-dependence
            for iFace in 1:my_mesh_2D.cellNodesCount[cellID]
                faceID = cell_faces[iFace]
                neighbor_cellID = cellNeighbors[iFace]

                #get the neighbor cell ID
                neighborID = my_mesh_2D.cellNeighbors_Dict[cellID][iFace]

                if !my_mesh_2D.bFace_is_boundary[faceID]  # If it is not a boundary face
                    #jac_sparsity[cellID, neighborID] = 1.0  # h does not directly depend on h of the neighbor cell
                    jac_sparsity[cellID+numOfCells, neighborID+numOfCells] = 1.0  # hu -> hu (neighbor)
                    jac_sparsity[cellID+2*numOfCells, neighborID+2*numOfCells] = 1.0  # hv -> hv (neighbor)
                end
            end          
        end
    end

    return jac_sparsity
end
