
#Problem constants and solution fields for 2D SWE with finite volume method

struct swe_2D_consts
    unit_system::String   #The choice of unit system, e.g., SI, EN
    k_n::Float64   #unit conversion factor for Manning's equation
    g::Float64   #gravity constant 
    t::Float64   #current time 
    dt_min::Float64  #minimum time step size 
    dt::Float64   #time step size 
    CFL::Float64   #CFL number 
    tStart::Float64  #start time for simulation
    tEnd::Float64  #end time for simulation 
    tspan::Tuple{Float64, Float64} #time span for simulation
    h_small::Float64  #dry bed water depth threshold (e.g., 1.0e-3)
    RiemannSolver::String   #The choice of Riemann solver, e.g., Roe, HLL, and HLLC 
end

