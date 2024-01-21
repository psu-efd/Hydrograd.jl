#semin-discretize the 1D SWE to convert it to an ODE 
#This file should be problem specific because we may invert different parameters 

# Problem-specific function to calculate the RHS of ODE: dQdt (= dhdt and dqdt)
# In this case, we will invert p = [inletQ, exitH, ManningN].
function swe_1D_rhs!(dQdt, Q, p, t, mesh, swe_1d_constants, left_bcType, right_bcType, S0)

    #for convenience
    dx = mesh.dx
    nCells = mesh.nCells
    nFaces = mesh.nFaces

    g = swe_1d_constants.g
    h_small = swe_1d_constants.h_small
    RiemannSolver = swe_1d_constants.RiemannSolver


    #fluxes on faces 
    fluxes = zeros(eltype(Q), 2, nFaces)

    #println("time t = ", t)

    #set the parameter values 
    left_bcValue = p[1] #hard-coded here. Should be more elegant. 
    right_bcValue = p[2]
    ManningN = p[3]

    #get a view (no need for a copy) of the current solution variables h and q 
    h = Q[:,1]
    q = Q[:,2]

    #ghost cells on the left and right boundaries
    #    wall (no flux, dh/dx=0 and q=hu=0), bcValue = 0.0 (not used)
    #    zeroGradient (dh/dx=0, dq/dx=0), bcValue = 0.0 (not used)
    #    inletQ (specify inlet specific discharge q=hu=q_in), bcValue = q_in
    #    exitH (specify outlet water depth, h = h_outlet), bcValue = h_outlet

    #left boundary
    if left_bcType=="wall"   #wall
        h_ghostCell_left = h[1]
        q_ghostCell_left = -q[1]    #flux goes the opposite direction -> no flux
    elseif left_bcType=="zeroGradient #zeroGradient"
        h_ghostCell_left = h[1]
        q_ghostCell_left = q[1]
    elseif left_bcType=="inletQ" #inletQ (bcValue=q_in)
        h_ghostCell_left = h[1]
        q_ghostCell_left = left_bcValue
    elseif left_bcType=="exitH" #exitH (bcValue=h_outlet)
        h_ghostCell_left = left_bcValue
        q_ghostCell_left = q[1]
    else 
        println("Left boundary condition type not recognized.")
        exit(-1)  #exit with an error code of -1
    end

    # right boundary
    if right_bcType == "wall"  # wall
        h_ghostCell_right = h[nCells]
        q_ghostCell_right = -q[nCells]
    elseif right_bcType == "zeroGradient"  # zeroGradient
        h_ghostCell_right = h[nCells]
        q_ghostCell_right = q[nCells]
    elseif right_bcType == "inletQ"  # inletQ (bcValue=q_in)
        h_ghostCell_right = h[nCells]
        q_ghostCell_right = right_bcValue
    elseif right_bcType == "exitH"  # exitH (bcValue=h_outlet)
        h_ghostCell_right = right_bcValue
        q_ghostCell_right = q[nCells]
    else 
        println("Right boundary condition type not recognized.")
        exit(-1)  #exit with an error code of -1
    end

    #loop through all faces
    @inbounds for iFace in 1:nFaces
        if (iFace == 1)               #left boundary face
            h_face = [h_ghostCell_left, h[1]]
            q_face = [q_ghostCell_left, q[1]]
            bcType = left_bcType
            bcValue = left_bcValue
        elseif (iFace == nFaces)   #right boundary face
            h_face = [h[nCells], h_ghostCell_right]
            q_face = [q[nCells], q_ghostCell_right]
            bcType = left_bcType
            bcValue = left_bcValue
        else                          #internal face
            h_face = [h[iFace-1], h[iFace]]
            q_face = [q[iFace-1], q[iFace]]
            bcType = "internal"
            bcValue = zero(eltype(p))  #not used
        end 

        if (RiemannSolver == "Roe")
            println("Not implemented yet")
            exit(-1)  #exit with an error code of -1
        elseif (RiemannSolver == "HLL")
            flux = Riemann_1D_hll(g, h_face, q_face, bcType, bcValue, h_small)
        elseif (RiemannSolver == "HLLC")
            println("Not implemented yet")
            exit(-1)  #exit with an error code of -1
        else
            println("Wrong choice of RiemannSolver")
            exit(-1)  #exit with an error code of -1
        end

        fluxes[1,iFace] = flux[1] 
        fluxes[2,iFace] = flux[2]
    end 

    #loop through all cells
    @inbounds for iCell in 1:nCells

        #calcuate the RHS of the ODEs
        #flux_east = @view fluxes[:,iCell+1] #flux on east face
        #flux_west = @view fluxes[:,iCell]   #flux on west face
        flux_east = fluxes[:,iCell+1] #flux on east face
        flux_west = fluxes[:,iCell]   #flux on west face

        #fields.dhdt[iCell] = - (flux_east[1]- flux_west[1]) / dx
        dQdt[iCell,1] = - (flux_east[1]- flux_west[1]) / dx

        if (h[iCell] <= h_small) #if a dry cell, no flow resistance term
            #fields.dqdt[iCell] = -(flux_east[2] - flux_west[2]) / dx + g * h[iCell] * fields.S0[iCell]
            dQdt[iCell,2] = -(flux_east[2] - flux_west[2]) / dx + g * h[iCell] * S0[iCell]
        else
            #fields.dqdt[iCell] = (- (flux_east[2] - flux_west[2]) / dx + g * h[iCell] * fields.S0[iCell]
            #           - g*ManningN^2/max(h[iCell], h_small)^(7.0/3.0)*abs(fields.q[iCell])*fields.q[iCell])
            dQdt[iCell,2] = (- (flux_east[2] - flux_west[2]) / dx + g * h[iCell] * S0[iCell]
                       - g*ManningN^2/max(h[iCell], h_small)^(7.0/3.0)*abs(q[iCell])*q[iCell])
        end
    end 

end
