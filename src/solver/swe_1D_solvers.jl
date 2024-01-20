#Solvers for 1D SWE 

#struct for Riemann state
#Ref: https://github.com/wme7/ShallowWaterEquations
Base.@kwdef mutable struct RiemannState
    L
    R
end

#struct for flux
#Ref: https://github.com/wme7/ShallowWaterEquations
Base.@kwdef mutable struct swe_1D_flux
    hFlux
    qFlux
end

# HLL Riemman solver
#Ref: https://github.com/wme7/ShallowWaterEquations
# bcType: -1 (for internal faces, not even a boundary face)
function Riemann_1D_hll(g, h_face, q_face, bcType, bcValue, h_small)

    h_L = h_face[1]
    h_R = h_face[2]
    q_L = q_face[1]
    q_R = q_face[2]

    if (h_L <= h_small && h_R <= h_small)
        if bcType==inletQ    #for inletQ BC, we need to specify the discharge (q=hu) even it is dry
            return [bcValue, zero(eltype(bcValue))]
        else
            return [zero(eltype(bcValue)), zero(eltype(bcValue))]
        end 
    end 

    if h_L <= h_small
        u_L = zero(eltype(bcValue))
    else
        u_L = q_L / h_L
    end

    if h_R <= h_small
        u_R = zero(eltype(bcValue))
    else
        u_R = q_R / h_R
    end

    if h_L <= h_small
        s_L = u_R - 2.0 * sqrt(g * h_R)
    else
        s_L = min(u_R - sqrt(g * h_R), u_L - sqrt(g * h_L))
    end 

    if h_R <= h_small
        s_R = u_L + 2.0 * sqrt(g * h_L)
    else
        s_R = max(u_R + sqrt(g * h_R), u_L + sqrt(g * h_L))
    end

    h_flux_L = q_L
    h_flux_R = q_R

    if (abs(s_R-s_L) < 1e-10)
        h_flux_star = zero(eltype(bcValue))
    else
        h_flux_star = (s_R * h_flux_L - s_L * h_flux_R + s_L * s_R * (h_R - h_L)) / (s_R - s_L)
    end 

    q_flux_L = q_L * u_L + g * h_L * h_L / 2.0
    q_flux_R = q_R * u_R + g * h_R * h_R / 2.0

    if (abs(s_R-s_L) < 1e-10)
        q_flux_star = zero(eltype(bcValue))
    else
        q_flux_star = (s_R * q_flux_L - s_L * q_flux_R + s_L * s_R * (q_R - q_L)) / (s_R - s_L)
    end 

    if (0.0 <= s_L)
        h_flux = h_flux_L
        q_flux = q_flux_L
    elseif (s_L <= 0 && 0 <= s_R)
        h_flux = h_flux_star
        q_flux = q_flux_star
    else
        h_flux = h_flux_R
        q_flux = q_flux_R
    end 

    return [h_flux, q_flux]
end 


# function to calculate the RHS of ODE: dQdt (= dhdt and dqdt)
function swe_1D_rhs!(dQdt, Q, p, t, mesh, fields, bcTypes, bcValues)

    #println("time t = ", t)

    #fluxex on faces 
    fluxes = zeros(eltype(Q), 2, mesh.nFaces)

    #unpack: p[1]=inletQ, p[2]=exitH, p[3]=ManningN, 
    #set the parameter values 

    left_bcType = bcTypes[1]
    right_bcType = bcTypes[2]
    left_bcValue = bcValues[1]
    right_bcValue = bcValues[2]

    #hard-coded here. Should be more elegant.
    left_bcValue = p[1]  
    right_bcValue = p[2]
    swe_1d_para.ManningN = p[3]
    #finished unpacking

    #get a view (no need for a copy) of the current solution variables h and q 
    h = Q[:,1]
    q = Q[:,2]

    #for convenience
    dx = mesh.dx
    g = fields.swe_1d_para.g 
    ManningN = fields.swe_1d_para.ManningN
    h_small = fields.swe_1d_para.h_small 
    RiemannSolver = fields.swe_1d_para.RiemannSolver

    nCells = mesh.nCells  #number of cells
    nFaces = mesh.nFaces  #number of faces

    #fluxes_on_faces = [] #np.zeros(2, N_faces) #python is row-major

    #ghost cells on the left and right boundaries
    #    wall (no flux, dh/dx=0 and q=hu=0), bcValue = 0.0 (not used)
    #    zeroGradient (dh/dx=0, dq/dx=0), bcValue = 0.0 (not used)
    #    inletQ (specify inlet specific discharge q=hu=q_in), bcValue = q_in
    #    exitH (specify outlet water depth, h = h_outlet), bcValue = h_outlet

    #left boundary
    if left_bcType=="wall"   #wall
        h_ghostCell_left = h[1]
        q_ghostCell_left = -q[1]    #flux goes the opposite direction -> no flux
    elseif left_bcType=="zeroGradient" #zeroGradient
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
    if right_bcType=="wall"
        h_ghostCell_right = h[nCells]
        q_ghostCell_right = -q[nCells]
    elseif right_bcType=="zeroGradient"  # zeroGradient
        h_ghostCell_right = h[nCells]
        q_ghostCell_right = q[nCells]
    elseif right_bcType=="inletQ"  # inletQ (bcValue=q_in)
        h_ghostCell_right = h[nCells]
        q_ghostCell_right = right_bcValue
    elseif right_bcType=="exitH"  # exitH (bcValue=h_outlet)
        h_ghostCell_right = right_bcValue
        q_ghostCell_right = q[nCells]
    else 
        println("Right boundary condition type not recognized.")
        exit(-1)  #exit with an error code of -1
    end

    #loop through all faces
    @inbounds for iFace in 1:mesh.nFaces
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
    for iCell in 1:mesh.nCells

        #calcuate the RHS of the ODEs
        flux_east = fluxes[:,iCell+1] #flux on east face
        flux_west = fluxes[:,iCell]   #flux on west face

        dQdt[iCell,1] = - (flux_east[1]- flux_west[1]) / dx

        if (h[iCell] <= h_small) #if a dry cell, no flow resistance term
            dQdt[iCell,2] = -(flux_east[2] - flux_west[2]) / dx + g * h[iCell] * S0[iCell]
        else
            dQdt[iCell,2] = (- (flux_east[2] - flux_west[2]) / dx + g * h[iCell] * fields.S0[iCell]
                       - g*ManningN^2/max(h[iCell], h_small)^(7.0/3.0)*abs(q[iCell])*q[iCell])
        end
    end 

end