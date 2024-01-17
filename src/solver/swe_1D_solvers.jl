#Solvers for 1D SWE 

#struct for Riemann state
#Ref: https://github.com/wme7/ShallowWaterEquations
Base.@kwdef mutable struct RiemannState
    L::Float64 = 0.0
    R::Float64 = 0.0
end

#struct for flux
#Ref: https://github.com/wme7/ShallowWaterEquations
Base.@kwdef mutable struct swe_1D_flux
    hFlux::Float64 = 0.0
    qFlux::Float64 = 0.0
end

# HLL Riemman solver
#Ref: https://github.com/wme7/ShallowWaterEquations
# bcType: -1 (for internal faces, not even a boundary face)
function Riemann_1D_hll(g, h_face, q_face, bcType, bcValue, h_small)
    h_L = h_face.L
    h_R = h_face.R
    q_L = q_face.L
    q_R = q_face.R

    if (h_L <= h_small && h_R <= h_small)
        if bcType==inletQ    #for inletQ BC, we need to specify the discharge (q=hu) even it is dry
            return swe_1D_flux(bcValue, 0.0)
        else
            return swe_1D_flux(0.0, 0.0)
        end 
    end 

    if h_L <= h_small
        u_L = 0.0
    else
        u_L = q_L / h_L
    end

    if h_R <= h_small
        u_R = 0.0
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
        h_flux_star = 0.0
    else
        h_flux_star = (s_R * h_flux_L - s_L * h_flux_R + s_L * s_R * (h_R - h_L)) / (s_R - s_L)
    end 

    q_flux_L = q_L * u_L + g * h_L * h_L / 2.0
    q_flux_R = q_R * u_R + g * h_R * h_R / 2.0

    if (abs(s_R-s_L) < 1e-10)
        q_flux_star = 0.0
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

    return swe_1D_flux(h_flux, q_flux)
end 


# function to calculate the RHS of ODE: dhdt and dqdt
function swe_1D_rhs!(mesh::mesh_1D, fields::swe_1D_fields)

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
    h_ghostCell_left = 0.0
    h_ghostCell_right = 0.0
    q_ghostCell_left = 0.0
    q_ghostCell_right = 0.0

    #left boundary
    if fields.leftBoundary.bcType==wall   #wall
        h_ghostCell_left = fields.h[1]
        q_ghostCell_left = -fields.q[1]    #flux goes the opposite direction -> no flux
    elseif fields.leftBoundary.bcType==zeroGradient #zeroGradient
        h_ghostCell_left = fields.h[1]
        q_ghostCell_left = fields.q[1]
    elseif fields.leftBoundary.bcType==inletQ #inletQ (bcValue=q_in)
        h_ghostCell_left = fields.h[1]
        q_ghostCell_left = fields.leftBoundary.bcValue
    elseif fields.leftBoundary.bcType==exitH #exitH (bcValue=h_outlet)
        h_ghostCell_left = fields.leftBoundary.bcValue
        q_ghostCell_left = fields.q[1]
    end

    # right boundary
    if fields.rightBoundary.bcType == wall  # wall
        h_ghostCell_right = fields.h[nCells]
        q_ghostCell_right = -fields.q[nCells]
    elseif fields.rightBoundary.bcType == zeroGradient  # zeroGradient
        h_ghostCell_right = fields.h[nCells]
        q_ghostCell_right = fields.q[nCells]
    elseif fields.rightBoundary.bcType == inletQ  # inletQ (bcValue=q_in)
        h_ghostCell_right = fields.h[nCells]
        q_ghostCell_right = fields.rightBoundary.bcValue
    elseif fields.rightBoundary.bcType == exitH  # exitH (bcValue=h_outlet)
        h_ghostCell_right = fields.rightBoundary.bcValue
        q_ghostCell_right = fields.q[nCells]
    end

    #loop through all faces
    for iFace in 1:nFaces
        if (iFace == 1)               #left boundary face
            h_face = RiemannState(h_ghostCell_left, fields.h[1])
            q_face = RiemannState(q_ghostCell_left, fields.q[1])
            bcType = fields.leftBoundary.bcType
            bcValue = fields.leftBoundary.bcValue
        elseif (iFace == nFaces)   #right boundary face
            h_face = RiemannState(fields.h[nCells], h_ghostCell_right)
            q_face = RiemannState(fields.q[nCells], q_ghostCell_right)
            bcType = fields.leftBoundary.bcType
            bcValue = fields.leftBoundary.bcValue
        else                          #internal face
            h_face = RiemannState(fields.h[iFace-1], fields.h[iFace])
            q_face = RiemannState(fields.q[iFace-1], fields.q[iFace])
            bcType = internal
            bcValue = 0.0  #not used
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

        fields.fluxes[1,iFace] = flux.hFlux 
        fields.fluxes[2,iFace] = flux.qFlux
    end 

    #loop through all cells
    for iCell in 1:nCells

        #calcuate the RHS of the ODEs
        flux_east = @view fields.fluxes[:,iCell+1] #flux on east face
        flux_west = @view fields.fluxes[:,iCell]   #flux on west face

        fields.dhdt[iCell] = - (flux_east[1]- flux_west[1]) / dx

        if (fields.h[iCell] <= h_small) #if a dry cell, no flow resistance term
            fields.dqdt[iCell] = -(flux_east[2] - flux_west[2]) / dx + g * fields.h[iCell] * fields.S0[iCell]
        else
            fields.dqdt[iCell] = (- (flux_east[2] - flux_west[2]) / dx + g * fields.h[iCell] * fields.S0[iCell]
                       - g*ManningN^2/max(fields.h[iCell], h_small)^(7.0/3.0)*abs(fields.q[iCell])*fields.q[iCell])
        end
    end 
end