

# using LinearAlgebra

# # Physical Constants
# const g = 9.81  # Gravitational acceleration

# # Flux Computation Using Roe or HLL Solver
# function compute_flux(hL, uL, vL, hR, uR, vR)
#     # Dry/wet conditions handling
#     if hL <= 1e-6 && hR <= 1e-6
#         return 0.0, 0.0, 0.0  # No flux if both sides are dry
#     end

#     # Roe or HLL solver (simplified example, extend as needed)
#     if hL > 1e-6 && hR > 1e-6  # Both wet
#         u_avg = 0.5 * (uL + uR)
#         h_avg = 0.5 * (hL + hR)
#         flux_h = h_avg * u_avg
#         flux_hu = flux_h * u_avg + 0.5 * g * h_avg^2
#         flux_hv = flux_h * vL  # Assumes aligned flow for simplicity
#         return flux_h, flux_hu, flux_hv
#     else  # Wet/dry interface
#         return 0.0, 0.0, 0.0  # Treat as a wall (adjust if needed)
#     end
# end

# # Apply Boundary Conditions
# function apply_boundary_conditions(U, boundary_faces, boundary_types)
#     for (face_id, btype) in boundary_faces
#         if btype == :inlet
#             U[face_id, :] .= [1.0, 0.1, 0.0]  # Example inlet flow discharge
#         elseif btype == :outlet
#             U[face_id, 1] = 1.0  # Fixed water surface elevation
#         elseif btype == :wall
#             U[face_id, 2:3] .= 0.0  # No-slip wall
#         end
#     end
# end

# # Finite Volume Update
# function update_cells(U, mesh, dt)
#     new_U = copy(U)
#     for cell_id in 1:size(mesh.cells, 1)
#         for face_id in mesh.cell_faces[cell_id]
#             hL, huL, hvL = U[cell_id, :]
#             hR, huR, hvR = U[mesh.neighbors[face_id], :]

#             # Convert to velocities
#             uL, vL = huL / hL, hvL / hL
#             uR, vR = huR / hR, hvR / hR

#             # Compute flux
#             flux_h, flux_hu, flux_hv = compute_flux(hL, uL, vL, hR, uR, vR)

#             # Update cell values
#             new_U[cell_id, 1] -= dt / mesh.cell_areas[cell_id] * flux_h
#             new_U[cell_id, 2] -= dt / mesh.cell_areas[cell_id] * flux_hu
#             new_U[cell_id, 3] -= dt / mesh.cell_areas[cell_id] * flux_hv
#         end
#     end
#     return new_U
# end

# # Main Solver Function
# function solve_shallow_water(mesh, U, boundary_faces, boundary_types, t_end, dt)
#     t = 0.0
#     while t < t_end
#         apply_boundary_conditions(U, boundary_faces, boundary_types)
#         U = update_cells(U, mesh, dt)
#         t += dt
#     end
#     return U
# end


# # Example mesh (replace with actual unstructured mesh data)
# mesh = (
#     cells = [1, 2, 3],                     # Example cell IDs
#     cell_faces = [[1, 2], [2, 3]],         # Faces for each cell
#     neighbors = [2, 1, 3],                 # Neighbors for each face
#     cell_areas = [1.0, 1.0, 1.0]           # Area of each cell
# )

# # Initial state (replace with realistic values)
# U = [
#     [1.0, 0.0, 0.0];  # Cell 1
#     [1.0, 0.0, 0.0];  # Cell 2
#     [1.0, 0.0, 0.0]   # Cell 3
# ]

# # Boundary conditions
# boundary_faces = Dict(1 => :inlet, 2 => :outlet, 3 => :wall)
# boundary_types = [:inlet, :outlet, :wall]

# # Solve
# t_end = 1.0
# dt = 0.01
# U_final = solve_shallow_water(mesh, U, boundary_faces, boundary_types, t_end, dt)

# println(U_final)


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

function roe_riemann_solver(hL, huL, hvL, hR, huR, hvR, g, n; hmin=1e-6)
    # Extract unit normal components
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

    # Compute Roe averages
    sqrt_hL = sqrt(hL)
    sqrt_hR = sqrt(hR)
    hRoe = sqrt_hL * sqrt_hR
    uRoe = (sqrt_hL * uL + sqrt_hR * uR) / (sqrt_hL + sqrt_hR)
    vRoe = (sqrt_hL * vL + sqrt_hR * vR) / (sqrt_hL + sqrt_hR)
    unRoe = uRoe * nx + vRoe * ny
    cRoe = sqrt(g * hRoe)

    # Compute wave speeds
    SL = unRoe - cRoe
    SR = unRoe + cRoe

    # Compute fluxes
    h_flux_L = huL * nx + hvL * ny
    hu_flux_L = (huL * uL + 0.5 * g * hL^2) * nx + huL * vL * ny
    hv_flux_L = (hvL * uL) * nx + (hvL * vL + 0.5 * g * hL^2) * ny

    h_flux_R = huR * nx + hvR * ny
    hu_flux_R = (huR * uR + 0.5 * g * hR^2) * nx + huR * vR * ny
    hv_flux_R = (hvR * uR) * nx + (hvR * vR + 0.5 * g * hR^2) * ny

    # Compute the Roe flux
    if SL >= 0
        # Use left state flux
        return h_flux_L, hu_flux_L, hv_flux_L
    elseif SR <= 0
        # Use right state flux
        return h_flux_R, hu_flux_R, hv_flux_R
    else
        # Roe flux
        h_flux = (SR * h_flux_L - SL * h_flux_R + SL * SR * (hR - hL)) / (SR - SL)
        hu_flux = (SR * hu_flux_L - SL * hu_flux_R + SL * SR * (huR - huL)) / (SR - SL)
        hv_flux = (SR * hv_flux_L - SL * hv_flux_R + SL * SR * (hvR - hvL)) / (SR - SL)
        return h_flux, hu_flux, hv_flux
    end
end

function roe_riemann_solver_with_zbed(hL, huL, hvL, hR, huR, hvR, zL, zR, g, n; hmin=1e-6)
    # Extract unit normal components
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

    # Compute Roe averages
    sqrt_hL = sqrt(hL)
    sqrt_hR = sqrt(hR)
    hRoe = (hL + hR)/2.0  #sqrt_hL * sqrt_hR
    uRoe = (sqrt_hL * uL + sqrt_hR * uR) / (sqrt_hL + sqrt_hR)
    vRoe = (sqrt_hL * vL + sqrt_hR * vR) / (sqrt_hL + sqrt_hR)
    unRoe = uRoe * nx + vRoe * ny
    cRoe = sqrt(g * hRoe)

    # Compute wave speeds
    SL = unRoe - cRoe
    SR = unRoe + cRoe

    # Bed slope term
    delta_z = zR - zL
    Sb_x = -g * hRoe * delta_z * nx
    Sb_y = -g * hRoe * delta_z * ny

    # Compute fluxes
    h_flux_L = huL * nx + hvL * ny
    hu_flux_L = (huL * uL + 0.5 * g * hL^2) * nx + huL * vL * ny
    hv_flux_L = (hvL * uL) * nx + (hvL * vL + 0.5 * g * hL^2) * ny

    h_flux_R = huR * nx + hvR * ny
    hu_flux_R = (huR * uR + 0.5 * g * hR^2) * nx + huR * vR * ny
    hv_flux_R = (hvR * uR) * nx + (hvR * vR + 0.5 * g * hR^2) * ny

    # Compute the Roe flux
    if SL >= 0
        # Use left state flux
        return h_flux_L, hu_flux_L + Sb_x, hv_flux_L + Sb_y
    elseif SR <= 0
        # Use right state flux
        return h_flux_R, hu_flux_R + Sb_x, hv_flux_R + Sb_y
    else
        # Roe flux
        h_flux = (SR * h_flux_L - SL * h_flux_R + SL * SR * (hR - hL)) / (SR - SL)
        hu_flux = (SR * hu_flux_L - SL * hu_flux_R + SL * SR * (huR - huL)) / (SR - SL) + Sb_x
        hv_flux = (SR * hv_flux_L - SL * hv_flux_R + SL * SR * (hvR - hvL)) / (SR - SL) + Sb_y
        return h_flux, hu_flux, hv_flux
    end
end


# Inputs
hL, huL, hvL = 2.0, 0.5, 0.0    # Left state: height and momentum
hR, huR, hvR = 1.0, 0.1, 0.0   # Right state: height and momentum
zL, zR = 0.0, 0.1               # Bed elevations (non-uniform bed)
g = 9.81                        # Gravitational constant
n = [1.0, 0.0]                  # Unit normal vector (face-aligned)

# Call the HLL Riemann solver
h_flux, hu_flux, hv_flux = hll_riemann_solver(hL, huL, hvL, hR, huR, hvR, g, n)

# Output results
println("hll solver:")
println("h_flux: $h_flux")
println("hu_flux: $hu_flux")
println("hv_flux: $hv_flux")

# Call the Roe Riemann solver
h_flux, hu_flux, hv_flux = roe_riemann_solver(hL, huL, hvL, hR, huR, hvR, g, n)

# Output results
println("roe solver:")
println("h_flux: $h_flux")
println("hu_flux: $hu_flux")
println("hv_flux: $hv_flux")

# Call the Roe Riemann solver with zbed
h_flux, hu_flux, hv_flux = roe_riemann_solver_with_zbed(hL, huL, hvL, hR, huR, hvR, zL, zR, g, n)

# Output results
println("roe solver zbed:")
println("h_flux: $h_flux")
println("hu_flux: $hu_flux")
println("hv_flux: $hv_flux")