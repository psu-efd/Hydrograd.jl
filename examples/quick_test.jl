using LinearAlgebra
using Printf

# Define constants
const g = 9.81  # Gravitational constant

# Define Riemann Solver
function riemann_solver(hL, hR, uL, uR, n)
    # Normal velocity components
    uLn = dot(uL, n)
    uRn = dot(uR, n)
    
    # Roe-averaged quantities
    h_avg = 0.5 * (hL + hR)
    u_avg = 0.5 * (uLn + uRn)
    c_avg = sqrt(g * h_avg)
    
    # Compute wave speeds
    sL = u_avg - c_avg
    sR = u_avg + c_avg
    
    # Flux splitting
    if sL >= 0
        return hL * uLn, hL * uLn^2 + 0.5 * g * hL^2
    elseif sR <= 0
        return hR * uRn, hR * uRn^2 + 0.5 * g * hR^2
    else
        h_flux = ((sR * hL - sL * hR) + sL * sR * (hR - hL)) / (sR - sL)
        momentum_flux = ((sR * (hL * uLn) - sL * (hR * uRn)) +
                         sL * sR * (hR * uRn - hL * uLn)) / (sR - sL)
        return h_flux, momentum_flux
    end
end

# Temporal integration using RK45
function rk45_step(mesh, h, u, Δt, compute_flux)
    # Define Butcher tableau for RK4 (standard)
    a = [0.0, 0.5, 0.5, 1.0]
    b = [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0]
    stages = length(a)
    
    # Stage calculations
    k_h = zeros(Float64, size(h, 1), stages)
    k_u = zeros(Float64, size(u, 1), stages, 2)
    h_temp = h
    u_temp = u
    
    for s in 1:stages
        k_h[:, s], k_u[:, s, :] = compute_flux(mesh, h_temp, u_temp)
        h_temp = h + Δt * a[s] * k_h[:, s]
        u_temp = u + Δt * a[s] * k_u[:, s, :]
    end
    
    # Final update
    h_new = h
    u_new = u
    for s in 1:stages
        h_new += Δt * b[s] * k_h[:, s]
        u_new += Δt * b[s] * k_u[:, s, :]
    end
    
    return h_new, u_new
end

# Compute flux for the entire mesh
function compute_flux(mesh, h, u)
    h_flux = zeros(size(h))
    u_flux = zeros(size(u))
    
    for face in mesh.faces
        cellL, cellR, n, area = face
        hL, hR = h[cellL], h[cellR]
        uL, uR = u[cellL, :], u[cellR, :]
        f_h, f_u = riemann_solver(hL, hR, uL, uR, n)
        h_flux[cellL] -= f_h * area
        h_flux[cellR] += f_h * area
        u_flux[cellL, :] -= f_u * area
        u_flux[cellR, :] += f_u * area
    end
    
    return h_flux, u_flux
end

# Main simulation loop
function simulate(mesh, h_init, u_init, t_final, Δt)
    h = h_init
    u = u_init
    t = 0.0
    
    while t < t_final
        h, u = rk45_step(mesh, h, u, Δt, compute_flux)
        t += Δt
        @printf("Time: %f\n", t)
    end
    
    return h, u
end

# Placeholder mesh data structure
mutable struct Mesh
    faces::Array{Tuple{Int, Int, Vector{Float64}, Float64}}
    nodes::Array{Vector{Float64}}
    cells::Array{Vector{Int}}
end

# Example usage
mesh = Mesh(
    faces = [
        (1, 2, [1.0, 0.0], 1.0), 
        (2, 3, [0.0, 1.0], 1.0)
    ],
    nodes = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
    cells = [[1, 2, 3]]
)

h_init = [1.0, 0.8, 0.9]
u_init = [0.0 0.0; 0.1 0.1; -0.1 -0.1]
t_final = 1.0
Δt = 0.01

h, u = simulate(mesh, h_init, u_init, t_final, Δt)

println(h)
println(u)

println("Simulation completed successfully")


# # Example usage
# # Define the mesh nodes
# nodes = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (2.0, 0.0), (2.0, 1.0)]

# # Define the cells (list of node indices)
# cells = [[1, 2, 3, 4], [2, 5, 6, 3]]

# # Define scalar variables (list of arrays, one for each variable)
# scalar_data = [[1.0, 2.0], [0.1, 0.2]]  # Example: Two scalar fields
# scalar_names = ["Pressure", "Temperature"]

# # Define vector variables (list of arrays, one for each variable)
# vector_data = [[[1.0, 0.0], [0.0, 1.0]], [[-1.0, 1.0], [1.0, -1.0]]]  # Example: Two vector fields
# vector_names = ["Velocity", "Force"]

# # Export to a VTK file
# file_path = joinpath(@__DIR__, "unstructured_mesh.vtk" ) 
# export_to_vtk(file_path, nodes, cells, scalar_data, scalar_names, vector_data, vector_names)

#println("VTK file written successfully")
