#Riemann solvers for 2D SWE: In fact, the Riemann problem is 1D in the normal direction of the face. 

#Roe Riemann solver
function Riemann_2D_Roe(iCell::Int, settings::ControlSettings, xiL::T1, hstillL::T2, hL::T3, huL::T4, hvL::T5, zb_L::T6, 
    xiR::T7, hstillR::T8, hR::T9, huR::T10, hvR::T11, zb_R::T12, zb_face::T13, 
    S0_on_face::Vector{T14},
    g::Float64, normal::Vector{T15}; hmin::Float64=1e-6) where {T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15}
    
    #Data type 
    data_type = promote_type(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15)

    # Extract unit normal components
    nx, ny = normal

    # Handle dry bed conditions (need to be smoothed; see notes)
    if (hL <= hmin && hR <= hmin)  #both sides are dry
        # Both sides are dry
        if settings.bVerbose
            print("Both sides are dry")
        end

        return zeros(data_type, 3)
    elseif (((hL + zb_L) < (zb_R + hmin)) && (hR <= hmin))    #left side WSE is below right side bed and right side is dry
        #In this case, the interface is virtually a wall. We will treat it as a wall.
        if settings.bVerbose
            print("Left side WSE is below right side bed and right side is dry")
        end

        #In this case, the interface is virtually a wall. We will treat it as a wall.
        #hR <- hL 
        #huR <- -huL
        #hvR <- -hvL

        hR = hL
        huR = -huL
        hvR = -hvL

        #return zeros(data_type, 3)
    elseif (((hR + zb_R) < (zb_L + hmin)) && (hL <= hmin))    #right side WSE is below left side bed and left side is dry
        if settings.bVerbose
            print("Right side WSE is below left side bed and left side is dry")
        end

        #In this case, the interface is virtually a wall. We will treat it as a wall.
        #hL <- hR 
        #huL <- -huR
        #hvL <- -hvR

        hL = hR
        huL = -huR
        hvL = -hvR

        #return zeros(data_type, 3)    
    elseif hL <= hmin                        #left side is dry, but right side is wet and WSE on right side is above left side bed
        # Left side is dry
        if settings.bVerbose
            print("Left side is dry")
        end
        
        h_flux = huR * nx + hvR * ny
        hu_flux = (huR * (huR / hR) + 0.5 * g * smooth_pow(hR, 2)) * nx + huR * (hvR / hR) * ny
        hv_flux = (hvR * (huR / hR)) * nx + (hvR * (hvR / hR) + 0.5 * g * smooth_pow(hR, 2)) * ny

        return [h_flux, hu_flux, hv_flux]
    elseif hR <= hmin                        #right side is dry, but left side is wet and WSE on left side is above right side bed
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
    R_mat = @SMatrix [data_type(0.0) data_type(1.0) data_type(1.0); 
                      ny uRoe-cRoe*nx uRoe+cRoe*nx; 
                      -nx vRoe-cRoe*ny vRoe+cRoe*ny]
    
    L_mat = @SMatrix [-(uRoe*ny-vRoe*nx) ny -nx; 
                      unRoe*over_two_cRoe+data_type(0.5) -nx*over_two_cRoe -ny*over_two_cRoe; 
                      -unRoe*over_two_cRoe+data_type(0.5) nx*over_two_cRoe ny*over_two_cRoe]
    
    absLamda = @SMatrix [smooth_abs(unRoe) data_type(0.0) data_type(0.0); 
                         data_type(0.0) smooth_abs(unRoe-cRoe) data_type(0.0); 
                         data_type(0.0) data_type(0.0) smooth_abs(unRoe+cRoe)]

    dQ = @SVector [xiR - xiL, huR - huL, hvR - hvL]

    #parenthesis to enforce only matrix-vector multiplication, not matrix-matrix multiplication, for better performance
    absA_dQ = R_mat * (absLamda * (L_mat * dQ))

    # Compute fluxes
    xi_flux_L = huL * nx + hvL * ny
    hu_flux_L = (huL * uL + 0.5 * g * (smooth_pow(xiL, 2) + 2.0*xiL*hstillL)) * nx + huL * vL * ny
    hv_flux_L = (hvL * uL) * nx + (hvL * vL + 0.5 * g * (smooth_pow(xiL, 2) + 2.0*xiL*hstillL)) * ny

    xi_flux_R = huR * nx + hvR * ny
    hu_flux_R = (huR * uR + 0.5 * g * (smooth_pow(xiR, 2) + 2.0*xiR*hstillR)) * nx + huR * vR * ny
    hv_flux_R = (hvR * uR) * nx + (hvR * vR + 0.5 * g * (smooth_pow(xiR, 2) + 2.0*xiR*hstillR)) * ny

    # Roe flux
    #absA_dQ = absA_dQ .* zero(data_type)  #for debugging
    xi_flux = (xi_flux_L + xi_flux_R - absA_dQ[1]) / 2.0   
    hu_flux = (hu_flux_L + hu_flux_R - absA_dQ[2]) / 2.0
    hv_flux = (hv_flux_L + hv_flux_R - absA_dQ[3]) / 2.0

    Zygote.ignore() do
        if iCell == -8  
            @show iCell
            println("before bed slope term:")
            @show xi_flux_L, hu_flux_L, hv_flux_L
            @show xi_flux_R, hu_flux_R, hv_flux_R
            @show absA_dQ
            @show xi_flux, hu_flux, hv_flux
            @show zb_L, zb_R
        end
    end


    #add the bed slope term contribution to the momentum fluxes    
    #slope_term = g * hRoe * (zb_L - zb_face)
    #hu_flux -= slope_term * nx
    #hv_flux -= slope_term * ny

    Zygote.ignore() do
        if iCell == -8  
            println("after bed slope term:")
            @show S0_on_face, hRoe
            @show slope_term
            @show xi_flux, hu_flux, hv_flux
        end
    end

    return [xi_flux, hu_flux, hv_flux]

end