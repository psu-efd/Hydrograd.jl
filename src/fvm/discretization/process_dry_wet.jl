#process dry wet flags
function process_dry_wet_flags(my_mesh_2D, h, zb_cells, swe_2D_constants)

    #cell wet/dry flag (0: dry, 1: wet)
    b_dry_wet = BitVector(h .> swe_2D_constants.h_small)

    # Process adjacent flags for all cells at once using array comprehension
    adjacent_flags = BitVector([
        let neighbors = my_mesh_2D.cellNeighbors_Dict[iCell]
            any(1:length(neighbors)) do iNeighbor
                neighbor = neighbors[iNeighbor]
                face_is_boundary = my_mesh_2D.bFace_is_boundary[my_mesh_2D.cellFacesList[iCell, iNeighbor]]
                
                face_is_boundary || 
                !b_dry_wet[neighbor]
            end
        end
        for iCell in 1:my_mesh_2D.numOfCells
    ])

    high_dry_flags = BitVector([
        let neighbors = my_mesh_2D.cellNeighbors_Dict[iCell]
            any(1:length(neighbors)) do iNeighbor
                neighbor = neighbors[iNeighbor]
                face_is_boundary = my_mesh_2D.bFace_is_boundary[my_mesh_2D.cellFacesList[iCell, iNeighbor]]
                
                face_is_boundary || 
                (!b_dry_wet[neighbor] && (h[iCell] + zb_cells[iCell]) < zb_cells[neighbor])
            end
        end
        for iCell in 1:my_mesh_2D.numOfCells
    ])

    return b_dry_wet, adjacent_flags, high_dry_flags
end