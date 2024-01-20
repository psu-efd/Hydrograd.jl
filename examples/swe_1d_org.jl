"""
This is a simple example for solving 1D SWE using Godnove scheme. 

This is the original code with my own ODE predictor-corrector solver. No SciML involed. 

The development of this 1D referenced the following Python implementation:
Ref: https://github.com/wme7/ShallowWaterEquations
"""

using Plots
using CSV
using DataFrames

using AdHydraulics

#println(@__FILE__)
#println(@__DIR__)

#Functions to setup the bed profile (this is just for some examples)
function my_gauss(x::Float64; sigma::Float64=1.0, h::Float64=1.0, mid::Float64=0.0)
    variance = sigma^2
    return h * exp(-(x-mid)^2/(2*variance))
end

function setup_bed!(mesh::mesh_1D, fields::swe_1D_fields)
    # parameters for bed setup
    b_bump_height = 0.3
    
    #loop through faces (points)
    for i in 1:mesh.nFaces
        fields.zb_face[i] = my_gauss(mesh.xFaces[i], sigma=1.0, h=b_bump_height, mid=mesh.L/2.0)
    end

    #loop through cells
    for i in 1:mesh.nCells
        fields.zb_cell[i] = (fields.zb_face[i+1] + fields.zb_face[i])/2.0
        fields.S0[i] = - (fields.zb_face[i+1] - fields.zb_face[i])/(mesh.xFaces[i+1] - mesh.xFaces[i])    
    end
end

# setup initial condition: free surface 
function setup_initial_eta!(mesh, fields)
    # parameters for initial free surface profile setup
    bump_center_x = 5.0  # center of the bump
    bump_half_width_x = 1.0  # bump's half width

    #loop over cells
    for i in 1:mesh.nCells
        #if abs(mesh.xCells[i] - bump_center_x) < bump_half_width_x:
        #    fields.eta[i] = 1.1
        #else:
        #    fields.eta[i] = 1.0

        if mesh.xCells[i] < bump_center_x
            fields.eta[i] = 1.0
        else
            fields.eta[i] = 0.5 #0.0  #0.5
        end
    end

    #update water depth
    fields.h = fields.eta - fields.zb_cell
    fields.h[fields.h .< 0.0] .= fields.swe_1d_para.h_small  #ensure positivity of water depth h

    #update the free surface elevation again in case h has been clipped
    fields.eta = fields.h + fields.zb_cell
end

#calculate the total water volume in the domain
function calc_total_water_volume(h, dx)
    total_water_volume = sum(h)*dx
    return total_water_volume
end 

#make plots of the results
function make_plots()
    data_dataframe = CSV.read(joinpath(dirname(@__FILE__), "results.csv"), DataFrame)
    data = Matrix(data_dataframe)  #x, h, q, eta, zb

    p1 = plot(data[:,1], data[:,4], label="free surface")
    plot!(data[:,1], data[:,5], label="bed")
    plot!(legend=:outerbottom, fg_legend = :transparent)

    data_dataframe = CSV.read(joinpath(dirname(@__FILE__), "total_water_volume.csv"), DataFrame)
    data = Matrix(data_dataframe)  #total_water_volume

    p2 = plot(data[:,1])
    
    display(p1)
    #display(p2)
end 

#create the 1D FVM mesh 
L = 10.0 #length of the domain
nCells = 200 #number of cells

mesh = initialize_mesh_1D(nCells, L)

#println(mesh.xCells)

#create 1D FVM fields 
#  create parameters
swe_1d_para = swe_1D_parameters(
    g = 9.81,   
    ManningN = 0.03,
    t = 0.0,   
    dt_min = 0.005, 
    dt = 0.005,
    CFL = 0.001,
    tEnd = 5.0,
    h_small = 1.0e-3,
    RiemannSolver = "HLL"
)

#  create the left and right boundaries
#leftBoundary = Boundary_1D(bcType=inletQ::Boundary_Type_Name, bcValue=0.5)
#rigthBoundary = Boundary_1D(bcType=exitH::Boundary_Type_Name, bcValue=1.0)

leftBoundary = Boundary_1D(bcType=wall::Boundary_Type_Name, bcValue=0.0)
rigthBoundary = Boundary_1D(bcType=wall::Boundary_Type_Name, bcValue=0.0)

fields = initialize_swe_1D_fields(mesh, swe_1d_para, leftBoundary, rigthBoundary)

# setup the bed 
setup_bed!(mesh, fields)

#plot the bed profile for checking 
#p1 = plot(mesh.xCells, fields.zb_cell, label="bed at cells")
#plot!(mesh.xFaces, fields.zb_face, label="bed at faces")
#p2 = plot(mesh.xCells, fields.S0, label="slope at cells")
#display(p1)
#display(p2)

#setup the initial free surface 
setup_initial_eta!(mesh, fields)

#plot the free surface profile for checking 
#p1 = plot(mesh.xCells, fields.zb_cell, label="bed at cells")
#plot!(mesh.xCells, fields.eta, label="free surface")
#plot!(mesh.xCells, fields.h, label="water depth")
#display(p1)

p = [0.5, 1.0, 0.03, mesh, fields]   #Solution parameter: inletQ, exitH, and ManningN 

#my own ODE solve 
while fields.swe_1d_para.t < fields.swe_1d_para.tEnd
    println("dt, t = ", fields.swe_1d_para.dt, fields.swe_1d_para.t)

    Q = hcat(fields.h, fields.q)

    # first Runge-Kutta (predictor)
    dQdt = swe_1D_rhs!(Q, p, fields.swe_1d_para.t)
    dhdt = dQdt[:,1]
    dqdt = dQdt[:,2]

    # solve old values
    h_old = copy(fields.h)
    q_old = copy(fields.q)
    dhdt_old = copy(dhdt)
    dqdt_old = copy(dqdt)
    for i in 1:mesh.nCells
        fields.h[i] = fields.h[i] + dhdt[i] * fields.swe_1d_para.dt
        fields.q[i] = fields.q[i] + dqdt[i] * fields.swe_1d_para.dt
    end 

    #clip_solution(N, h, q, h_small)

    Q = hcat(fields.h, fields.q)

    # second Runge-Kutta (corrector)
    dQdt = swe_1D_rhs!(Q, p, fields.swe_1d_para.t)
    dhdt = dQdt[:,1]
    dqdt = dQdt[:,2]

    for i in 1:mesh.nCells
        fields.h[i] = h_old[i] + (dhdt[i] + dhdt_old[i]) * fields.swe_1d_para.dt / 2.0
        fields.q[i] = q_old[i] + (dqdt[i] + dqdt_old[i]) * fields.swe_1d_para.dt / 2.0
    end 

    #clip_solution(N, h, q, h_small)

    Q = hcat(fields.h, fields.q)

    #update free surface elevation eta 
    fields.eta = fields.h + fields.zb_cell

    #calculate CFL to determine time step size
    #println("dt_cfl = ", fields.swe_1d_para.CFL * mesh.dx / maximum(abs.(fields.q)./fields.h))
    #fields.swe_1d_para.dt = max(fields.swe_1d_para.dt_min, fields.swe_1d_para.CFL * 
    #                             mesh.dx / maximum(abs.(fields.q)./fields.h))

    push!(fields.total_water_volume, calc_total_water_volume(fields.h, mesh.dx))

    fields.swe_1d_para.t = fields.swe_1d_para.t + fields.swe_1d_para.dt
end 

#save results to files the same directory of this jl file
open(joinpath(dirname(@__FILE__), "results.csv"), "w") do fo
    println(fo, "x, h, q, eta, zb")
    for i in 1:mesh.nCells
      println(fo, mesh.xCells[i], ", ", fields.h[i], ", ", fields.q[i], ", ", fields.eta[i], ", ", fields.zb_cell[i])
    end
  end

open(joinpath(dirname(@__FILE__), "total_water_volume.csv"), "w") do fo
    println(fo, "total_water_volume")
    for volume in fields.total_water_volume
      println(fo, volume)
    end
  end

make_plots()

println("Done!")
