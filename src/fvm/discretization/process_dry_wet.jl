#process dry wet flags
function process_dry_wet_flags(my_mesh_2D, h, zb_cells, swe_2D_constants)

    #cell wet/dry flag (0: dry, 1: wet)
    b_dry_wet = h .> swe_2D_constants.h_small  

    #initialize the flags for whether a cell is adjacent to dry land and high dry land
    b_Adjacent_to_dry_land = falses(my_mesh_2D.numOfCells)
    b_Adjacent_to_high_dry_land = falses(my_mesh_2D.numOfCells)

    #loop over all cells: check each of the neighbors of the current cell, if any of the neighbors is dry or wall boundary, the current cell is adjacent to dry land
    for iCell in 1:my_mesh_2D.numOfCells

        #get the neighbors of the current cell
        neighbors = my_mesh_2D.cellNeighbors_Dict[iCell]

        #@show iCell, neighbors

        #loop over all the neighbors
        for iNeighbor in eachindex(neighbors)

            #if the neighbor is a wall boundary, the current cell is adjacent to dry land. Both flags are set to true
            if my_mesh_2D.bFace_is_boundary[my_mesh_2D.cellFacesList[iCell, iNeighbor]]
                b_Adjacent_to_dry_land[iCell] = true
                b_Adjacent_to_high_dry_land[iCell] = true
            else   #if the neighbor is not a wall boundary, check if the neighbor is a dry cell
                if b_dry_wet[neighbors[iNeighbor]] == 0.0
                    b_Adjacent_to_dry_land[iCell] = true
                end

                #check if the neighbor is high dry land: WSE of current cell is lower than the zb of the neighbor
                if ~b_dry_wet[neighbors[iNeighbor]] && (h[iCell] + zb_cells[iCell]) < zb_cells[neighbors[iNeighbor]]
                    b_Adjacent_to_high_dry_land[iCell] = true
                end
            end
        end

    end

    return b_dry_wet, b_Adjacent_to_dry_land, b_Adjacent_to_high_dry_land
end