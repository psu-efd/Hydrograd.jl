# process bed profiles, such as setup bed profile, interplation of bed elevation from face to cell,
# cell to face, calculate slope S0, etc. 
#
# This file should be problem specific (each problem may have different bed profile)

#Functions to setup the bed profile (this is just for some examples)
function my_gauss(x::Float64; sigma::Float64=1.0, h::Float64=1.0, mid::Float64=0.0)
    variance = sigma^2
    return h * exp(-(x - mid)^2 / (2 * variance))
end

function setup_bed!(mesh, zb_face, zb_cell, S0)

    # parameters for bed setup
    b_bump_height = 0.3

    #loop through cells
    @inbounds for i in 1:mesh.nCells
        zb_cell[i] = my_gauss(mesh.xCells[i], sigma=1.0, h=b_bump_height, 
                              mid=maximum(mesh.xFaces) / 2.0)
    end

    interploate_zb_from_cell_to_face_and_compute_S0!(mesh, zb_face, zb_cell, S0)

    # optionally plot the bed for checking
    #plot_bed(mesh, zb_cell, zb_face, S0)
end

function interploate_zb_from_cell_to_face_and_compute_S0!(mesh, zb_face, zb_cell, S0)

    #loop through faces to compute zb_face 
    @inbounds for i in 1:mesh.nFaces
        if i==1    #left boundary face 
            zb_face[i] = zb_cell[i]
        elseif i==mesh.nFaces   #right boundary face 
            zb_face[i] = zb_cell[i-1]
        else
            zb_face[i] = (zb_cell[i-1] + zb_cell[i]) / 2.0
        end 
    end

    #loop through cells to compute S0 at cells 
    @inbounds for i in 1:mesh.nCells
        S0[i] = -(zb_face[i+1] - zb_face[i]) / (mesh.xFaces[i+1] - mesh.xFaces[i])
    end
end

function interploate_zb_from_face_to_cell_and_compute_S0!(mesh, zb_face, zb_cell, S0)
    #loop through cells
    @inbounds for i in 1:mesh.nCells
        zb_cell[i] = (zb_face[i+1] + zb_face[i]) / 2.0
        S0[i] = -(zb_face[i+1] - zb_face[i]) / (mesh.xFaces[i+1] - mesh.xFaces[i])
    end
end

function interploate_zb_from_face_to_cell!(mesh, zb_face, zb_cell_local)

    #loop through cells
    @inbounds for i in 1:mesh.nCells
        zb_cell_local[i] = (zb_face[i+1] + zb_face[i]) / 2.0
    end

end

#plot the bed profile for checking 
function plot_bed(mesh, zb_cell, zb_face, S0)
    p1 = plot(mesh.xCells, zb_cell, label="bed at cells")
    plot!(mesh.xFaces, zb_face, label="bed at faces")
    p2 = plot(mesh.xCells, S0, label="slope at cells")
    display(p1)
    display(p2)
end
