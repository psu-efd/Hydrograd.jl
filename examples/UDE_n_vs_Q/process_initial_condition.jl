# process initial condition, such as setting up initial free surface, water depth, etc. 
# This file should be problem specific because each problem should have different ICs. 

# setup initial condition: free surface 
function setup_initial_eta!(mesh, eta, zb_cell, h, swe_1d_constants)
    # parameters for initial free surface profile setup
    bump_center_x = 5.0  # center of the bump
    bump_half_width_x = 1.0  # bump's half width

    nCells = mesh.nCells
    xCells = mesh.xCells
    h_small = swe_1d_constants.h_small

    #loop over cells
    @inbounds for i in 1:nCells
        #if abs(mesh.xCells[i] - bump_center_x) < bump_half_width_x:
        #    fields.eta[i] = 1.1
        #else:
        #    fields.eta[i] = 1.0

        if xCells[i] < bump_center_x
            eta[i] = 1.0
        else
            eta[i] = 0.5 #0.5  #0.5
        end
    end

    #update water depth
    @inbounds for i in 1:nCells
        h[i] = eta[i] - zb_cell[i]
    end

    h[h.<0.0] .= h_small  #ensure positivity of water depth h

    #update the free surface elevation again in case h has been clipped
    @inbounds for i in 1:nCells
        eta[i] = h[i] + zb_cell[i]
    end

    #optionally plot the free surface for checking 
    #plot_free_surface_elevation(xCells, zb_cell, eta, h)
end

#plot the free surface profile for checking 
function plot_free_surface_elevation(xCells, zb_cell, eta, h)
    p1 = plot(xCells, zb_cell, label="bed at cells")
    plot!(xCells, eta, label="free surface")
    plot!(xCells, h, label="water depth")
    display(p1)
end
