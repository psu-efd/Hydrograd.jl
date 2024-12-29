
#Problem constants and solution fields for 1D SWE with finite volume method

Base.@kwdef mutable struct swe_1D_consts
    g::Float64 = 9.81   #gravity constant 

    t::Float64   #time 
    dt_min::Float64 = 0.001  #minimum time step size 
    dt::Float64   #time step size 
    CFL::Float64 = 0.4  #CFL number 
    tEnd::Float64  #end time for simulation 

    h_small::Float64 = 1.0e-3  #dry bed water depth threshold (e.g., 1.0e-3)

    RiemannSolver::String = "HLL"  #The choice of Riemann solver, e.g., Roe, HLL, and HLLC 
end