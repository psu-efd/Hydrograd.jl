#Riemann solvers for 2D SWE: In fact, the Riemann problem is 1D in the normal direction of the face. 

#Roe Riemann solver
function Riemann_2D_Roe(settings::ControlSettings, hL::T1, huL::T2, hvL::T3, hR::T4, huR::T5, hvR::T6, 
    g::Float64, normal::Vector{T7}; hmin::Float64=1e-6) where {T1, T2, T3, T4, T5, T6, T7}
    
    #Data type 
    data_type = promote_type(T1, T2, T3, T4, T5, T6, T7)

    # Extract unit normal components
    nx, ny = normal

    # Handle dry bed conditions (need to be smoothed; see notes)
    if hL < hmin && hR < hmin
        # Both sides are dry
        if settings.bVerbose
            print("Both sides are dry")
        end

        return zeros(data_type, 3)
    elseif hL < hmin
        # Left side is dry
        if settings.bVerbose
            print("Left side is dry")
        end
        
        h_flux = huR * nx + hvR * ny
        hu_flux = (huR * (huR / hR) + 0.5 * g * smooth_pow(hR, 2)) * nx + huR * (hvR / hR) * ny
        hv_flux = (hvR * (huR / hR)) * nx + (hvR * (hvR / hR) + 0.5 * g * smooth_pow(hR, 2)) * ny

        return [h_flux, hu_flux, hv_flux]
    elseif hR < hmin
        # Right side is dry
        if settings.bVerbose
            print("Right side is dry")
        end
        
        h_flux = huL * nx + hvL * ny
        hu_flux = (huL * (huL / hL) + 0.5 * g * smooth_pow(hL, 2)) * nx + huL * (hvL / hL) * ny
        hv_flux = (hvL * (huL / hL)) * nx + (hvL * (hvL / hL) + 0.5 * g * smooth_pow(hL, 2)) * ny

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
    sqrt_hL = smooth_sqrt(hL)
    sqrt_hR = smooth_sqrt(hR)
    hRoe = (hL + hR)/2.0    #Arithmetic average
    uRoe = (sqrt_hL * uL + sqrt_hR * uR) / (sqrt_hL + sqrt_hR)
    vRoe = (sqrt_hL * vL + sqrt_hR * vR) / (sqrt_hL + sqrt_hR)
    unRoe = uRoe * nx + vRoe * ny
    cRoe = smooth_sqrt(g * hRoe)
    over_two_cRoe = one(data_type)/2.0/cRoe

    # Compute wave speeds
    SL = unRoe - cRoe
    SR = unRoe + cRoe

    # define matrices 
    R_mat = @SMatrix [data_type(0.0) data_type(1.0) data_type(1.0); ny uRoe-cRoe*nx uRoe+cRoe*nx; -nx vRoe-cRoe*ny vRoe+cRoe*ny]
    #L_mat = [-(uRoe*ny-vRoe*nx) ny -nx; unRoe/2/cRoe+0.5 -nx/2/cRoe -ny/2/cRoe; -unRoe/2/cRoe+0.5 nx/2/cRoe ny/2/cRoe]
    L_mat = @SMatrix [-(uRoe*ny-vRoe*nx) ny -nx; unRoe*over_two_cRoe+data_type(0.5) -nx*over_two_cRoe -ny*over_two_cRoe; -unRoe*over_two_cRoe+data_type(0.5) nx*over_two_cRoe ny*over_two_cRoe]
    absLamda = @SMatrix [smooth_abs(unRoe) data_type(0.0) data_type(0.0); data_type(0.0) smooth_abs(unRoe-cRoe) data_type(0.0); data_type(0.0) data_type(0.0) smooth_abs(unRoe+cRoe)]

    #absA = R_mat * absLamda * L_mat 

    #dQ = @SVector (hR-hL, huR-huL, hvR-hvL)
    dQ = @SVector [hR - hL, huR - huL, hvR - hvL]

    #absA_dQ = absA * dQ

    #parenthesis to enforce only matrix-vector multiplication, not matrix-matrix multiplication, for better performance
    absA_dQ = R_mat * (absLamda * (L_mat * dQ))

    # Compute fluxes
    h_flux_L = huL * nx + hvL * ny
    hu_flux_L = (huL * uL + 0.5 * g * smooth_pow(hL, 2)) * nx + huL * vL * ny
    hv_flux_L = (hvL * uL) * nx + (hvL * vL + 0.5 * g * smooth_pow(hL, 2)) * ny

    h_flux_R = huR * nx + hvR * ny
    hu_flux_R = (huR * uR + 0.5 * g * smooth_pow(hR, 2)) * nx + huR * vR * ny
    hv_flux_R = (hvR * uR) * nx + (hvR * vR + 0.5 * g * smooth_pow(hR, 2)) * ny

    # Roe flux
    h_flux = (h_flux_L + h_flux_R - absA_dQ[1]) / 2.0   
    hu_flux = (hu_flux_L + hu_flux_R - absA_dQ[2]) / 2.0
    hv_flux = (hv_flux_L + hv_flux_R - absA_dQ[3]) / 2.0

    return [h_flux, hu_flux, hv_flux]

end