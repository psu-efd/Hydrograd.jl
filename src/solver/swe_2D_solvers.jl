#Riemann solvers for 2D SWE: In fact, the Riemann problem is 1D in the normal direction of the face. 

# HLL Riemman solver
#Ref: https://github.com/wme7/ShallowWaterEquations
# bcType: -1 (for internal faces, not even a boundary face)
function Riemann_2D_hll!(flux, g, h_face, q_face, bcType, bcValue, h_small)

    h_L = h_face[1]
    h_R = h_face[2]
    q_L = q_face[1]
    q_R = q_face[2]

    if (h_L <= h_small && h_R <= h_small)
        if bcType=="inletQ"    #for inletQ BC, we need to specify the discharge (q=hu) even it is dry
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
        if (h_L < 0.0 || h_R < 0.0)
            println("h_L, h_R = ", h_L, ", ", h_R)
        end
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

    flux[1] = h_flux
    flux[2] = q_flux

    #return [h_flux, q_flux]
end 

function hll_riemann_solver(hL, huL, hvL, hR, huR, hvR, g, n; hmin=1e-6)
    # Extract normal components
    nx, ny = n

    # Handle dry bed conditions
    if hL < hmin && hR < hmin
        # Both sides are dry
        return 0.0, 0.0, 0.0
    elseif hL < hmin
        # Left side is dry
        h_flux = huR * nx + hvR * ny
        hu_flux = (huR * (huR / hR) + 0.5 * g * hR^2) * nx + huR * (hvR / hR) * ny
        hv_flux = (hvR * (huR / hR)) * nx + (hvR * (hvR / hR) + 0.5 * g * hR^2) * ny
        return h_flux, hu_flux, hv_flux
    elseif hR < hmin
        # Right side is dry
        h_flux = huL * nx + hvL * ny
        hu_flux = (huL * (huL / hL) + 0.5 * g * hL^2) * nx + huL * (hvL / hL) * ny
        hv_flux = (hvL * (huL / hL)) * nx + (hvL * (hvL / hL) + 0.5 * g * hL^2) * ny
        return h_flux, hu_flux, hv_flux
    end

    # Compute velocities on left and right
    uL = huL / hL
    vL = hvL / hL
    uR = huR / hR
    vR = hvR / hR

    # Compute normal velocities
    unL = uL * nx + vL * ny
    unR = uR * nx + vR * ny

    # Compute wave speeds (assuming shallow water equations)
    cL = sqrt(g * hL)
    cR = sqrt(g * hR)

    SL = min(unL - cL, unR - cR)
    SR = max(unL + cL, unR + cR)

    # Handle cases for SL and SR
    if SL >= 0
        # Use left state flux
        h_flux = huL * nx + hvL * ny
        hu_flux = (huL * uL + 0.5 * g * hL^2) * nx + huL * vL * ny
        hv_flux = (hvL * uL) * nx + (hvL * vL + 0.5 * g * hL^2) * ny
        return h_flux, hu_flux, hv_flux
    elseif SR <= 0
        # Use right state flux
        h_flux = huR * nx + hvR * ny
        hu_flux = (huR * uR + 0.5 * g * hR^2) * nx + huR * vR * ny
        hv_flux = (hvR * uR) * nx + (hvR * vR + 0.5 * g * hR^2) * ny
        return h_flux, hu_flux, hv_flux
    else
        # Use HLL flux
        h_flux_L = huL * nx + hvL * ny
        h_flux_R = huR * nx + hvR * ny
        hu_flux_L = (huL * uL + 0.5 * g * hL^2) * nx + huL * vL * ny
        hu_flux_R = (huR * uR + 0.5 * g * hR^2) * nx + huR * vR * ny
        hv_flux_L = (hvL * uL) * nx + (hvL * vL + 0.5 * g * hL^2) * ny
        hv_flux_R = (hvR * uR) * nx + (hvR * vR + 0.5 * g * hR^2) * ny

        h_flux = (SR * h_flux_L - SL * h_flux_R + SL * SR * (hR - hL)) / (SR - SL)
        hu_flux = (SR * hu_flux_L - SL * hu_flux_R + SL * SR * (huR - huL)) / (SR - SL)
        hv_flux = (SR * hv_flux_L - SL * hv_flux_R + SL * SR * (hvR - hvL)) / (SR - SL)
        return h_flux, hu_flux, hv_flux
    end
end

function Riemann_2D_Roe(hL, huL, hvL, hR, huR, hvR, g, normal; hmin=1e-6)
    
    #Data type 
    T = eltype(hL)

    # Extract unit normal components
    nx, ny = normal

    # Handle dry bed conditions
    if hL < hmin && hR < hmin
        # Both sides are dry
        print("Both sides are dry")
        
        return zeros(eltype(hL), 3)
    elseif hL < hmin
        # Left side is dry
        print("Left side is dry")
        h_flux = huR * nx + hvR * ny
        hu_flux = (huR * (huR / hR) + 0.5 * g * hR^2) * nx + huR * (hvR / hR) * ny
        hv_flux = (hvR * (huR / hR)) * nx + (hvR * (hvR / hR) + 0.5 * g * hR^2) * ny

        return [h_flux, hu_flux, hv_flux]
    elseif hR < hmin
        # Right side is dry
        print("Right side is dry")
        h_flux = huL * nx + hvL * ny
        hu_flux = (huL * (huL / hL) + 0.5 * g * hL^2) * nx + huL * (hvL / hL) * ny
        hv_flux = (hvL * (huL / hL)) * nx + (hvL * (hvL / hL) + 0.5 * g * hL^2) * ny

        return [h_flux, hu_flux, hv_flux]
    end

    # Compute velocities on left and right
    uL = huL / hL
    vL = hvL / hL
    uR = huR / hR
    vR = hvR / hR

    # Compute normal velocities
    unL = uL * nx + vL * ny
    unR = uR * nx + vR * ny

    # Compute Roe averages
    sqrt_hL = sqrt(hL)
    sqrt_hR = sqrt(hR)
    hRoe = (hL + hR)/2.0    #Arithmetic average
    uRoe = (sqrt_hL * uL + sqrt_hR * uR) / (sqrt_hL + sqrt_hR)
    vRoe = (sqrt_hL * vL + sqrt_hR * vR) / (sqrt_hL + sqrt_hR)
    unRoe = uRoe * nx + vRoe * ny
    cRoe = sqrt(g * hRoe)

    # Compute wave speeds
    SL = unRoe - cRoe
    SR = unRoe + cRoe

    # define matrices 
    R_mat = [T(0.0) T(1.0) T(1.0); ny uRoe-cRoe*nx uRoe+cRoe*nx; -nx vRoe-cRoe*ny vRoe+cRoe*ny]
    L_mat = [-(uRoe*ny-vRoe*nx) ny -nx; unRoe/2/cRoe+0.5 -nx/2/cRoe -ny/2/cRoe; -unRoe/2/cRoe+0.5 nx/2/cRoe ny/2/cRoe]
    absLamda = [abs(unRoe) T(0.0) T(0.0); T(0.0) abs(unRoe-cRoe) T(0.0); T(0.0) T(0.0) abs(unRoe+cRoe)]

    absA = R_mat * absLamda * L_mat 

    dQ = [hR-hL; huR-huL; hvR-hvL]

    absA_dQ = absA * dQ

    # Compute fluxes
    h_flux_L = huL * nx + hvL * ny
    hu_flux_L = (huL * uL + 0.5 * g * hL^2) * nx + huL * vL * ny
    hv_flux_L = (hvL * uL) * nx + (hvL * vL + 0.5 * g * hL^2) * ny

    h_flux_R = huR * nx + hvR * ny
    hu_flux_R = (huR * uR + 0.5 * g * hR^2) * nx + huR * vR * ny
    hv_flux_R = (hvR * uR) * nx + (hvR * vR + 0.5 * g * hR^2) * ny

    # Roe flux
    h_flux = (h_flux_L + h_flux_R - absA_dQ[1]) / 2.0   
    hu_flux = (hu_flux_L + hu_flux_R - absA_dQ[2]) / 2.0
    hv_flux = (hv_flux_L + hv_flux_R - absA_dQ[3]) / 2.0

    return [h_flux, hu_flux, hv_flux]

end