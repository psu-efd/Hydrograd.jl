#Riemann solvers for 1D SWE 

# HLL Riemman solver
#Ref: https://github.com/wme7/ShallowWaterEquations
# bcType: -1 (for internal faces, not even a boundary face)
function Riemann_1D_hll!(flux, g, h_face, q_face, bcType, bcValue, h_small)

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