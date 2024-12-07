using Revise

using JLD2

using IterTools

using PyCall

using AdHydraulics

using OrdinaryDiffEq

#SciML
using SciMLSensitivity

#Optimizers
using Optimization
using OptimizationOptimisers
#using OptimizationPolyalgorithms
#using OptimizationNLopt
#using Optim
#using OptimizationFlux
#using Flux

#AD engines
#using Zygote
#using ForwardDiff
#using ReverseDiff

#for Bayesian estimation
#using Turing

# Load StatsPlots for visualizations and diagnostics.
#using StatsPlots

#using LinearAlgebra

#using Random
#Random.seed!(1234)

include("process_SRH_2D_input.jl")
include("process_ManningN.jl")
include("preprocess_BCs.jl")
include("process_bed_2D.jl")
include("process_initial_condition_2D.jl")
include("semi_discretize_2D.jl")
include("misc_utilities_2D.jl")

println("Solving 2D SWE...")

#define control variables
bSimulate_Synthetic_Data = false    #whether to do the 1D SWE simulation to create synthetic data 
bPlot_Simulation_Results = false   #whether to plot simulation results
bPerform_Inversion = true           #whether to do inversion 
bPlot_Inversion_Results = false     #whehter to plot the inversion results

#options for inversion 
bInversion_slope_loss = true   #whether to include slope loss 
bInversion_include_u = true    #whehter to add velocity to the loss function (by default, we already have water surface elevatino eta)

#directory to save to (the same directory as the main file)
save_path = dirname(@__FILE__)

#define a swe_2D_const
swe_2D_constants = swe_2D_consts(t=0.0, dt=0.1, tStart=0.0, tEnd=1.0)

#read data from SRH-2D hydro, geom, and material files
#srhhydro_file_name = "simple.srhhydro"
srhhydro_file_name = "oneD_channel_with_bump.srhhydro"
#srhhydro_file_name = "twoD_channel_with_bump.srhhydro"

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

#setup bed elevation: computer zb at cell centers from nodes, then interpolate zb from cell to face and compute the bed slope at cells
setup_bed!(my_mesh_2D.numOfCells, my_mesh_2D.numOfNodes, my_mesh_2D.nodeCoordinates, 
my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, my_mesh_2D.cell_centroids, zb_cells, zb_ghostCells, zb_faces, S0, true)

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

#preprocess boundary conditions 
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
    inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, 
    inletQ_ManningN, inletQ_Length, inletQ_TotalA, inletQ_DryWet)

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

#println("inletQ_faceIDs: ", inletQ_faceIDs)
#println("exitH_faceIDs: ", exitH_faceIDs)
#println("wall_faceIDs: ", wall_faceIDs)
#println("symm_faceIDs: ", symm_faceIDs)

#set up ODE parameters. In this problem, the "parameter" is para=zb_cells.
para = zb_cells
Q0 = hcat(h, q_x, q_y)   #initial condition: Q = [h q_x q_y]

Q_ghost = hcat(h_ghostCells, q_x_ghostCells, q_y_ghostCells)   #ghost cell values of Q

# time information 
#tspan = (swe_2D_constants.tStart, swe_2D_constants.tEnd)
#dt = swe_2D_constants.dt
#t = tspan[1]:dt:tspan[2]

tspan = (0.0, 10.0)    #100.0
dt = 0.05
t = tspan[1]:dt:tspan[2]

dt_save = (tspan[2] - tspan[1])/10.0
t_save = tspan[1]:dt_save:tspan[2]

# # define the ODE
ode_f = ODEFunction((dQdt, Q, para, t) -> swe_2D_rhs!(dQdt, Q, Q_ghost, para, t, my_mesh_2D, zb_cells, zb_ghostCells, zb_faces, S0,
            swe_2D_constants, ManningN_cells, ManningN_ghostCells,
            nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
            inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length,
            inletQ_TotalA, inletQ_DryWet,  
            nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
            exitH_internalCellIDs, exitH_faceCentroids, exitH_WSE,
            nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, 
            wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals,
            nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs, 
            symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals
            ),
            jac_prototype=nothing
        )

prob = ODEProblem(ode_f, Q0, tspan, para)

if bSimulate_Synthetic_Data

    solver_choice = "SciML"
    #solver_choice = "MyOwn"

    println("   Performing 2D SWE simulation ...")

    #for direct sensivity analysis
    #prob = ODEForwardSensitivityProblem((Q, p, t) -> swe_1D_rhs!(Q, p, t, 
    #                          mesh, swe_1d_constants, left_bcType, right_bcType, S0), Q0, tspan, p)

    # solve the ODE with a solver of your choice 
    #bm = @benchmark solve(prob, Heun(), adaptive=false, dt=dt, saveat=t)
    #@show bm 

    if solver_choice == "SciML"

        println("   Performing 2D SWE simulation with SciML solver ...")
    
        sol = solve(prob, Tsit5(), adaptive=false, dt=dt, saveat=t_save)
        #sol = solve(prob, Heun(), adaptive=false, dt=dt, saveat=t)
        #sol = solve(prob, AB3(), adaptive=false, dt=dt, saveat = t)
        #sol = solve(prob, Rosenbrock23(), adaptive=false, dt=dt, saveat = t)  #implicit

        #x, dp = extract_local_sensitivities(sol)
        #da = dp[1]
        #plot(sol.t, da', lw = 3)

        #println("Press any key to exit ...")
        #readline()
        #exit()

        # #save the results
        # #save the simulation solution results
        jldsave(joinpath(save_path, "simulation_solution.jld2"); sol)

        swe_2D_save_results(sol, total_water_volume, my_mesh_2D, zb_cells, save_path)
    end

    #My own ODE Solver
    if solver_choice == "MyOwn"

        println("   Performing 2D SWE simulation with MyOwn solver ...")

        # Finite Volume Update
        function update_cells!(dQdt, Q, my_mesh_2D, zb_cells, zb_ghostCells, zb_faces, S0, t, dt)
            
            # compute dQdt 
            swe_2D_rhs!(dQdt, Q, Q_ghost, para, t, my_mesh_2D, zb_cells, zb_ghostCells, zb_faces, S0,
            swe_2D_constants, ManningN_cells, ManningN_ghostCells,
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
        function solve_shallow_water(my_mesh_2D, dQdt, Q, zb_cells, zb_ghostCells, zb_faces, S0, t_end, dt)
            t = 0.0
            iStep = 0
            while t < t_end
                println("t = ", t)
                
                if iStep == 5000
                    print("Here.")
                end
                
                update_cells!(dQdt, Q, my_mesh_2D, zb_cells, zb_ghostCells, zb_faces, S0, t, dt)
                
                if iStep % 100 == 0
                    #compute and record the total water volume
                    push!(total_water_volume, [t, sum(Q[:,1] .* my_mesh_2D.cell_areas)])
                end
                
                if true && iStep % 200 == 0
                    
                    u_temp = Q[:,2] ./ Q[:,1]
                    v_temp = Q[:,3] ./ Q[:,1]
                    U_vector = hcat(u_temp, v_temp)
                    
                    vector_data = [U_vector] 
                    vector_names = ["U"]
                    
                    scalar_data = [Q[:,1], Q[:,2], Q[:,3], eta, zb_cells, ManningN_cells]
                    scalar_names = ["h", "hu", "hv", "WSE", "zb_cell", "ManningN"]
                    
                    file_path = joinpath(save_path, "solution_$(iStep).vtk" ) 
                    export_to_vtk_2D(file_path, my_mesh_2D.nodeCoordinates, my_mesh_2D.cellNodesList, my_mesh_2D.cellNodesCount, scalar_data, scalar_names, vector_data, vector_names)    
                end 
                
                t += dt
                iStep += 1
            end
            
        end
        
        t_end = 200.0
        dt = 0.05
        
        dQdt = zeros(Float64, my_mesh_2D.numOfCells, 3)
        
        solve_shallow_water(my_mesh_2D, dQdt, Q0, zb_cells, zb_ghostCells, zb_faces, S0, t_end, dt)
        
        #save total water volume to file 
        open(joinpath(save_path, "total_water_volume.csv"), "w") do fo
            println(fo, "time, total_water_volume")
            for time_volume in total_water_volume
                println(fo, time_volume[1], ",", time_volume[2])
            end
        end
        
    end

end

if bPlot_Simulation_Results

    println("   Plotting 2D SWE simulation results ...")

    #open the simulation result
    sol = load(joinpath(save_path, "simulation_solution.jld2"))["sol"]

    #plot the results at tEnd
    #swe_1D_make_plots(save_path)

    #make an animation of the simulation as a function of time 
    #swe_1D_make_animation(sol, mesh, zb_cell, save_path)
end


################
#Inversion part# 
################

if bPerform_Inversion

    println("   Performing inversion ...")

    #open the simulation result (as the ground truth)
    sol = load(joinpath(save_path, "simulation_solution.jld2"))["sol"]
    h_truth = Array(sol)[:,1,end]
    u_truth = Array(sol)[:,2,end]./Array(sol)[:,1,end]
    v_truth = Array(sol)[:,3,end]./Array(sol)[:,1,end]
    #@show sol

    #inversion parameters: zb_cell
    # initial guess for the parameters
    ps = zeros(my_mesh_2D.numOfCells)

    function predict(θ)
        #Array(solve(prob, Heun(), adaptive=false, p=θ, dt=dt, saveat=t))[:,1,end]
        Array(solve(prob, Tsit5(), adaptive=false, p=θ, dt=dt, saveat=t_save))  #[:,1,end]
    end

    SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true
    #jac = ForwardDiff.jacobian(predict, ps)
    #jac = Zygote.jacobian(predict, ps)
    #jac = ReverseDiff.jacobian(predict, ps)
    #@show jac
    #plot(jac)
    #exit()

    ## Defining Loss function
    #zb_cell_local = zeros(Number, mesh.nCells)

    function loss(θ)
        pred = predict(θ)            #Forward prediction with current θ (=zb at cells)
        l = pred[:,1,end] - h_truth  #loss = free surface elevation mismatch
        loss_pred_eta = sum(abs2, l)

        loss_pred_uv = 0.0

        if bInversion_include_u      #if also include u in the loss 
            l_u = pred[:,2,end]./pred[:,1,end] .- u_truth
            l_v = pred[:,3,end]./pred[:,1,end] .- v_truth

            loss_pred_uv = sum(abs2, l_u) + sum(abs2, l_v)
        end 

        loss_pred = loss_pred_eta + loss_pred_uv

        loss_slope = 0.0

        if bInversion_slope_loss    #if bed slope is included in the loss 
            loss_slope = calc_slope_loss(θ, my_mesh_2D)
        end 

        loss_total = loss_pred + loss_slope

        return loss_total, loss_pred, loss_pred_eta, loss_pred_uv, loss_slope, pred
    end

    #grad = Zygote.gradient(loss, ps)
    #grad = ForwardDiff.gradient(loss, ps)
    #grad = ForwardDiff.gradient(predict, θ)
    #@show grad
    #println("grad = ", grad)
    #readline()
    #exit()


    LOSS = []                              # Loss accumulator
    PRED = []                              # prediction accumulator
    PARS = []                              # parameters accumulator

    callback = function (θ, loss_total, loss_pred, loss_pred_eta, loss_pred_uv, loss_slope, pred) #callback function to observe training
        iter = size(LOSS)[1]  #get the inversion iteration number (=length of LOSS array)
        println("      iter, loss_total, loss_pred, loss_pred_eta, loss_pred_uv, loss_slope = ", iter, ", ", 
                  loss_total, ", ", loss_pred, ", ", loss_pred_eta, ", ", loss_pred_uv, ", ", loss_slope)

        append!(PRED, [pred[:,1,end]])
        append!(LOSS, [[loss_total, loss_pred, loss_pred_eta, loss_pred_uv, loss_slope, pred]])

        if !isa(θ, Vector{Float64})  #NLopt returns an optimization object, not an arrary
            #println("theta.u = ", θ.u)
            append!(PARS, [copy(θ.u)])
        else
            append!(PARS, [θ])
        end

        #if l > 1e-9
        #    false
        #else 
        #    true   #force the optimizer to stop 
        #end

        false
    end

    # Let see prediction vs. Truth
    #scatter(sol[:, end], label = "Truth", size = (800, 500))
    #plot!(PRED[end][:, end], lw = 2, label = "Prediction")

    #The following is from SciMLSensitivity documentation regarding the choice of AD 
    #  AutoForwardDiff(): The fastest choice for small optimizations
    #  AutoReverseDiff(compile=false): A fast choice for large scalar optimizations
    #  AutoTracker(): Like ReverseDiff but GPU-compatible
    #  AutoZygote(): The fastest choice for non-mutating array-based (BLAS) functions
    #  AutoFiniteDiff(): Finite differencing, not optimal but always applicable
    #  AutoModelingToolkit(): The fastest choice for large scalar optimizations
    #  AutoEnzyme(): Highly performant AD choice for type stable and optimized code

    adtype = Optimization.AutoZygote()         #works on MBP, but slow. Not on Windows (don't know why). 
    
    #adtype = Optimization.AutoForwardDiff()  #ForwardDiff.jl is not supported for this problem.
    #adtype = Optimization.AutoReverseDiff()
    #adtype = Optimization.AutoEnzyme()   #failed, don't use "not implemented yet error".

    optf = Optimization.OptimizationFunction((θ, p) -> loss(θ), adtype)

    #setup the bounds for the bathymetry
    lb_p = zeros(my_mesh_2D.numOfCells)
    lb_p .= -0.1  

    ub_p = zeros(my_mesh_2D.numOfCells)
    ub_p .= 0.3 

    #optprob = Optimization.OptimizationProblem(optf, ps, lb=lb_p, ub=ub_p)
    optprob = Optimization.OptimizationProblem(optf, ps)

    #res = Optimization.solve(optprob, PolyOpt(), callback = callback)  #PolyOpt does not support lb and ub 
    #res = Optimization.solve(optprob, NLopt.LD_LBFGS(), callback = callback)   #very fast 
    #res = Optimization.solve(optprob, Optim.BFGS(), callback=callback; iterations=30, maxiters=40, f_calls_limit=20, show_trace=true)
    #res = Optimization.solve(optprob, Optim.BFGS(), callback=callback, maxiters = 100; show_trace=false)  #f_tol=1e-3, iterations=10, local_maxiters = 10
    #res = Optimization.solve(optprob, Optim.LBFGS(), callback=callback)  #oscilates around 1e-7
    #res = Optimization.solve(optprob, Optim.Newton(), callback=callback)  #error: not supported as the Fminbox optimizer
    #res = Optimization.solve(optprob, Optim.GradientDescent(), callback=callback)  #very slow decrease in loss 
    res = Optimization.solve(optprob, Adam(0.01), callback=callback, maxiters=100)
    
    @show res
    @show res.original

    #save the inversion results
    #@show PARS 
    jldsave(joinpath(save_path, "inversion_results.jld2"); LOSS, PRED, PARS, sol)

end


println("All done!")
