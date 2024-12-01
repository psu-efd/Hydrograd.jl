
#Problem constants and solution fields for 2D SWE with finite volume method

Base.@kwdef mutable struct swe_2D_consts
    g::Float64 = 9.81   #gravity constant 

    t::Float64   #current time 
    dt_min::Float64 = 0.001  #minimum time step size 
    dt::Float64   #time step size 
    CFL::Float64 = 0.4  #CFL number 
    tStart::Float64  #start time for simulation
    tEnd::Float64  #end time for simulation 

    h_small::Float64 = 1.0e-3  #dry bed water depth threshold (e.g., 1.0e-3)

    RiemannSolver::String = "Roe"  #The choice of Riemann solver, e.g., Roe, HLL, and HLLC 
end