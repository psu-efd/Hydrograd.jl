# process bed profiles, such as setup bed profile, interplation of bed elevation from face to cell,
# cell to face, calculate slope S0, etc. 
#
# This file should be problem specific (each problem may have different bed profile)

#loss due to slope regularization
function calc_slope_loss(zb_cell, dx)
    """
    2nd order Central difference for 1st degree derivative
    """
    #[[zero(eltype(zb_cell))]; (zb_cell[3:end] - zb_cell[1:(end - 2)]) ./ (2.0 * dx); [zero(eltype(zb_cell))]]

    S0_max = 0.2
    S0_center = 0.0

    #return max.(S0_center, (slope_temp .- S0_max)) 
    return sum(max.(S0_center, abs.([(zb_cell[3:end] - zb_cell[1:(end - 2)]) ./ (2.0 * dx)][1]) .- S0_max))
end

#Functions to setup the bed elevation (this is just for some examples)
function my_gauss(x::Float64; sigma::Float64=1.0, h::Float64=1.0, mid::Float64=0.0)
    variance = sigma^2
    return h * exp(-(x - mid)^2 / (2 * variance))
end

#Set up the bed elevation
function setup_bed!(numOfCells, numOfNodes, nodeCoordinates, cellNodesList, cellNodesCount, cell_centroids, zb_cell, bPlotBed::Bool=false)
    # parameters for bed setup
    b_bump_height = 0.0

    xMid = (minimum(nodeCoordinates[:,1])+maximum(nodeCoordinates[:,1])) / 2.0

    #loop through cells
    @inbounds for i in 1:numOfCells
        zb_cell[i] = my_gauss(cell_centroids[i,1], sigma=1.0, h=b_bump_height, 
                              mid=xMid)
    end

    #interploate_zb_from_cell_to_face_and_compute_S0!(mesh, zb_face, zb_cell, S0)

    # optionally plot the bed for checking
    if bPlotBed
        vector_data = [] 
        vector_names = []

        scalar_data = [zb_cell]
        scalar_names = ["zb_cell"]

        file_path = joinpath(@__DIR__, "bed.vtk" ) 
        export_to_vtk_2D(file_path, nodeCoordinates, cellNodesList, cellNodesCount, scalar_data, scalar_names, vector_data, vector_names)    
        println("Bed is saved to ", file_path)
    end
    
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
