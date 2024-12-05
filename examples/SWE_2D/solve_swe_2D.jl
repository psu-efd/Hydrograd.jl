using Revise

using JLD2

using IterTools

using PyCall

using AdHydraulics

using OrdinaryDiffEq

include("process_SRH_2D_input.jl")
include("process_ManningN.jl")
include("preprocess_BCs.jl")
include("process_bed_2D.jl")
include("process_initial_condition_2D.jl")
include("semi_discretize_2D.jl")
include("misc_utilities_2D.jl")

println("Solving 2D SWE...")

#directory to save to (the same directory as the main file)
save_path = dirname(@__FILE__)

#define a swe_2D_const
swe_2D_constants = swe_2D_consts(t=0.0, dt=0.1, tStart=0.0, tEnd=1.0)

#read data from SRH-2D hydro, geom, and material files
srhhydro_file_name = "oneD_channel_with_bump.srhhydro"
srh_all_Dict = process_SRH_2D_input!(srhhydro_file_name)

#update swe_2D_constants based on the SRH-2D data
update_swe_2D_constants!(swe_2D_constants, srh_all_Dict)

#get the 2D mesh 
my_mesh_2D = srh_all_Dict["my_mesh_2D"]

#  setup bed elevation 
zb_faces = zeros(Float64, my_mesh_2D.numOfFaces)      #zb at faces 
zb_cells = zeros(Float64, my_mesh_2D.numOfCells)      #zb at cell centers 
zb_ghostCells = zeros(Float64, my_mesh_2D.numOfAllBounaryFaces)   #zb at ghost cell centers 
S0 = zeros(Float64, my_mesh_2D.numOfCells, 2)          #bed slope at cell centers 

setup_bed!(my_mesh_2D.numOfCells, my_mesh_2D.numOfNodes, my_mesh_2D.nodeCoordinates, 
my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, my_mesh_2D.cell_centroids, zb_cells, true)

#setup zb_ghostCells: zb_ghostCells[i] has the same bed elevation at the neighbor cell 
update_ghost_cells_scalar!(my_mesh_2D.numOfAllBounaryFaces, my_mesh_2D.allBoundaryFacesIDs_List, 
my_mesh_2D.faceCells_Dict, zb_cells, zb_ghostCells)

#interpolate zb from cell to face 
cells_to_faces_scalar!(my_mesh_2D.numOfFaces, my_mesh_2D.faceCells_Dict, zb_cells, zb_faces)

#compute bed slope at cell centers
compute_scalar_gradients!(my_mesh_2D.numOfCells, my_mesh_2D.cell_areas, my_mesh_2D.cell_normals, 
my_mesh_2D.face_lengths, my_mesh_2D.cellNodesCount, my_mesh_2D.cellFacesList, my_mesh_2D.cellNeighbors_Dict, 
zb_cells, S0)

#bed slope is negative of zb gradient 
S0 = -S0

println("zb_cells: ", zb_cells)   
println("zb_faces: ", zb_faces)
println("S0: ", S0) 

if true
    vector_data = [S0] 
    vector_names = ["S0"]
    
    scalar_data = [zb_cells]
    scalar_names = ["zb_cells"]
    
    file_path = joinpath(@__DIR__, "zb_S0.vtk" ) 
    export_to_vtk_2D(file_path, my_mesh_2D.nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, 
    scalar_data, scalar_names, vector_data, vector_names)    
    
    println("zb and S0 are saved to ", file_path)
    #exit(0)
end

#make a copy of the bathymetry truth: zb_cell_truth
zb_cells_truth = deepcopy(zb_cells)

#setup Manning's n 
ManningN_cells = zeros(my_mesh_2D.numOfCells)                    #Manning's n at cells 
ManningN_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)    #Manning's n at ghost cells

setup_ManningN!(ManningN_cells, ManningN_ghostCells, my_mesh_2D, srh_all_Dict)

# setup initial conditions 
eta = zeros(my_mesh_2D.numOfCells)          #free surface elevation at cells 
h = zeros(my_mesh_2D.numOfCells)            #water depth at cells 
q_x = zeros(my_mesh_2D.numOfCells)          #q_x=hu at cells 
q_y = zeros(my_mesh_2D.numOfCells)          #q_y=hv at cells 

eta_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #free surface elevation at ghost cells 
h_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)            #water depth at ghost cells 
q_x_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #q_x=hu at ghost cells 
q_y_ghostCells = zeros(my_mesh_2D.numOfAllBounaryFaces)          #q_y=hv at ghost cells 

total_water_volume = []   #total volume of water in the domain 

#setup initial condition for eta, h, q_x, q_y
setup_initial_eta!(my_mesh_2D.numOfCells, my_mesh_2D.nodeCoordinates, my_mesh_2D.cellNodesList, 
my_mesh_2D.cellNodesCount, my_mesh_2D.cell_centroids, eta, zb_cells, h, true)

#setup ghost cells for initial condition
update_ghost_cells_eta_h_q!(my_mesh_2D.numOfAllBounaryFaces, my_mesh_2D.allBoundaryFacesIDs_List, my_mesh_2D.faceCells_Dict, 
eta, h, q_x, q_y, eta_ghostCells, h_ghostCells, q_x_ghostCells, q_y_ghostCells)

#process boundary conditions 
inletQ_BC_indices = Int[]   #indices of inlet-Q boundaries in the boundary list
exitH_BC_indices = Int[]    #indices of exit-H boundaries in the boundary list
wall_BC_indices = Int[]     #indices of wall boundaries in the boundary list
symm_BC_indices = Int[]     #indices of symmetry (no slip) boundaries in the boundary list

#compute the indices of each boundary in the global boundary list (nWall_BCs and nSymm_BCs are updated in this function)
compute_boundary_indices!(my_mesh_2D, srh_all_Dict, inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices)

nInletQ_BCs = srh_all_Dict["nInletQ_BCs"]   #number of inlet-Q boundaries
nExitH_BCs = srh_all_Dict["nExitH_BCs"]     #number of exit-H boundaries
nWall_BCs = srh_all_Dict["nWall_BCs"]       #number of wall boundaries
nSymm_BCs = srh_all_Dict["nSymm_BCs"]       #number of symmetry boundaries

#some data arrays for inlet-q bounaries:
inletQ_faceIDs = Array{Array{Int}}(undef, nInletQ_BCs)   #face IDs of the inlet-q boundaries
inletQ_ghostCellIDs = Array{Array{Int}}(undef, nInletQ_BCs)   #ghost cell IDs of the inlet-q boundaries
inletQ_internalCellIDs = Array{Array{Int}}(undef, nInletQ_BCs)   #internal cell IDs of the inlet-q boundaries

inletQ_faceCentroids = Vector{Matrix{Float64}}(undef, nInletQ_BCs)   #face centroids of the inlet-q boundaries
inletQ_faceOutwardNormals = Vector{Matrix{Float64}}(undef, nInletQ_BCs)   #face outward normals of the inlet-q boundaries

inletQ_TotalQ = zeros(Float64, nInletQ_BCs)   #total discharge for each inlet-q boundary
inletQ_H = Array{Array{Float64}}(undef, nInletQ_BCs)   #inlet water depth for each face in each inlet-q boundary
inletQ_A = Array{Array{Float64}}(undef, nInletQ_BCs)   #area for each face in each inlet-q boundary
inletQ_ManningN = Array{Array{Float64}}(undef, nInletQ_BCs)   #Manning's n for each face in each inlet-q boundary
inletQ_Length = Array{Array{Float64}}(undef, nInletQ_BCs)   #length for each face in each inlet-q boundary
inletQ_TotalA = zeros(Float64, nInletQ_BCs)   #total cross-sectional area for each inlet-q boundary
inletQ_DryWet = Array{Array{Int}}(undef, nInletQ_BCs)   #dry(=0)/wet(=1) flag for each inlet-q boundary

#preprocess inlet-q boundaries
preprocess_inlet_q_boundaries(my_mesh_2D, nInletQ_BCs, srh_all_Dict["srhhydro_IQParams"], inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
    inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length, inletQ_TotalA, inletQ_DryWet)

#some data arrays for exit-h bounaries:
exitH_faceIDs = Array{Array{Int}}(undef, nExitH_BCs)   #face IDs of the exit-h boundaries
exitH_ghostCellIDs = Array{Array{Int}}(undef, nExitH_BCs)   #ghost cell IDs of the exit-h boundaries
exitH_internalCellIDs = Array{Array{Int}}(undef, nExitH_BCs)   #internal cell IDs of the exit-h boundaries

exitH_faceCentroids = Vector{Matrix{Float64}}(undef, nExitH_BCs)   #face centroids of the exit-h boundaries
exitH_faceOutwardNormals = Vector{Matrix{Float64}}(undef, nExitH_BCs)   #face outward normals of the exit-h boundaries

exitH_WSE = Array{Float64}(undef, nExitH_BCs)   #WSE for each exit-h boundary

exitH_H = Array{Array{Float64}}(undef, nExitH_BCs)   #inlet water depth for each exit-h boundary
exitH_A = Array{Array{Float64}}(undef, nExitH_BCs)   #inlet cross-sectional area for each exit-h boundary

#preprocess exit-h boundaries
preprocess_exit_h_boundaries(my_mesh_2D, nExitH_BCs, srh_all_Dict["srhhydro_EWSParamsC"], exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
    exitH_internalCellIDs, exitH_faceCentroids, exitH_faceOutwardNormals, exitH_WSE, exitH_H, exitH_A)

#some data arrays for wall bounaries:
wall_faceIDs = Array{Array{Int}}(undef, nWall_BCs)   #face IDs of the wall boundaries
wall_ghostCellIDs = Array{Array{Int}}(undef, nWall_BCs)   #ghost cell IDs of the wall boundaries
wall_internalCellIDs = Array{Array{Int}}(undef, nWall_BCs)   #internal cell IDs of the wall boundaries

wall_faceCentroids = Vector{Matrix{Float64}}(undef, nWall_BCs)   #face centroids of the wall boundaries
wall_outwardNormals = Vector{Matrix{Float64}}(undef, nWall_BCs)   #face outward normals of the wall boundaries

wall_H = Array{Array{Float64}}(undef, nWall_BCs)   #water depth for each wall boundary
wall_A = Array{Array{Float64}}(undef, nWall_BCs)   #cross-sectional area for each wall boundary

#preprocess wall boundaries
preprocess_wall_boundaries(my_mesh_2D, nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs,
    wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals, wall_H, wall_A)

#some data arrays for symmetry boundaries:
symm_faceIDs = Array{Array{Int}}(undef, nSymm_BCs)
symm_ghostCellIDs = Array{Array{Int}}(undef, nSymm_BCs)
symm_internalCellIDs = Array{Array{Int}}(undef, nSymm_BCs)

symm_faceCentroids = Vector{Matrix{Float64}}(undef, nSymm_BCs)
symm_outwardNormals = Vector{Matrix{Float64}}(undef, nSymm_BCs)

symm_H = Array{Array{Float64}}(undef, nSymm_BCs)
symm_A = Array{Array{Float64}}(undef, nSymm_BCs)

#preprocess symmetry boundaries
preprocess_symmetry_boundaries(my_mesh_2D, nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs,
    symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals, symm_H, symm_A)

println("inletQ_faceIDs: ", inletQ_faceIDs)
println("exitH_faceIDs: ", exitH_faceIDs)
println("wall_faceIDs: ", wall_faceIDs)
println("symm_faceIDs: ", symm_faceIDs)

#set up ODE parameters. In this problem, the "parameter" is para=zb_cells.
para = zb_cells
Q0 = hcat(h, q_x, q_y)   #initial condition: Q = [h q_x q_y]

Q_ghost = hcat(h_ghostCells, q_x_ghostCells, q_y_ghostCells)   #ghost cell values of Q

# time information 
#tspan = (swe_2D_constants.tStart, swe_2D_constants.tEnd)
#dt = swe_2D_constants.dt
#t = tspan[1]:dt:tspan[2]

tspan = (0.0, 1.0)
dt = 0.1
t = tspan[1]:dt:tspan[2]

dt_save = (tspan[2] - tspan[1])/10.0
t_save = tspan[1]:dt_save:tspan[2]

# # define the ODE
# ode_f = ODEFunction((dQdt, Q, para, t) -> swe_2D_rhs!(dQdt, Q, Q_ghost, para, t, my_mesh_2D, swe_2D_constants, ManningN_cells, ManningN_ghostCells,
#         nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
#         inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length,
#         inletQ_TotalA, inletQ_DryWet,  
#         nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
#         exitH_internalCellIDs, exitH_faceCentroids, exitH_WSE,
#         nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, 
#         wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals
#         ),
#     jac_prototype=nothing
#     )

# prob = ODEProblem(ode_f, Q0, tspan, para)

# sol = solve(prob, Tsit5(), adaptive=true, dt=dt, saveat=t)

# #save the results
# #save the simulation solution results
# jldsave(joinpath(save_path, "simulation_solution.jld2"); sol)

#swe_2D_save_results(sol, total_water_volume, my_mesh_2D, zb_cells, save_path)

# Finite Volume Update
function update_cells!(dQdt, Q, my_mesh_2D, t, dt)
    
    # compute dQdt 
    swe_2D_rhs!(dQdt, Q, Q_ghost, para, t, my_mesh_2D, swe_2D_constants, ManningN_cells, ManningN_ghostCells,
         nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
         inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length,
         inletQ_TotalA, inletQ_DryWet,  
         nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
         exitH_internalCellIDs, exitH_faceCentroids, exitH_WSE,
         nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, 
         wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals,
         nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs, 
         symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals
         )

    for iCell in 1:my_mesh_2D.numOfCells
        Q[iCell, :] += dt * dQdt[iCell, :]        

        #update water depth in case of dry cells
        if Q[iCell, 1] < swe_2D_constants.h_small
            Q[iCell, 1] = swe_2D_constants.h_small
            Q[iCell, 2] = 0.0
            Q[iCell, 3] = 0.0
        end

        #update WSE
        eta[iCell] = Q[iCell, 1] + zb_cells[iCell]
    end
    

end

# Main Solver Function
function solve_shallow_water(my_mesh_2D, dQdt, Q, t_end, dt)
    t = 0.0
    iStep = 0
    while t < t_end
        println("t = ", t)

        if iStep == 5000
            print("Here.")
        end

        update_cells!(dQdt, Q, my_mesh_2D, t, dt)

        if iStep % 100 == 0
           #compute and record the total water volume
            push!(total_water_volume, [t, sum(Q[:,1] .* my_mesh_2D.cell_areas)])
        end

        if true && iStep % 1000 == 0

            u_temp = Q[:,2] ./ Q[:,1]
            v_temp = Q[:,3] ./ Q[:,1]
            U_vector = hcat(u_temp, v_temp)

            vector_data = [U_vector] 
            vector_names = ["U"]
            
            scalar_data = [Q[:,1], Q[:,2], Q[:,3], eta, zb_cells, ManningN_cells]
            scalar_names = ["h", "hu", "hv", "eta", "zb_cell", "ManningN"]
            
            file_path = joinpath(save_path, "solution_$(iStep).vtk" ) 
            export_to_vtk_2D(file_path, my_mesh_2D.nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, scalar_data, scalar_names, vector_data, vector_names)    
        end 

        t += dt
        iStep += 1
    end
    return Q
end

t_end = 1000.0
dt = 0.01

dQdt = zeros(Float64, my_mesh_2D.numOfCells, 3)

solve_shallow_water(my_mesh_2D, dQdt, Q0, t_end, dt)

#save total water volume to file 
open(joinpath(save_path, "total_water_volume.csv"), "w") do fo
    println(fo, "time, total_water_volume")
    for time_volume in total_water_volume
        println(fo, time_volume[1], ",", time_volume[2])
    end
end


println("All done!")
