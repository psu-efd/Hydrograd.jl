
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

#update swe_2d_constants based on the SRH-2D data
function update_swe_2D_constants!(swe_2D_constants, srh_all_Dict)

    srhhydro_SimTime = srh_all_Dict["srhhydro_SimTime"]  
    tStart = srhhydro_SimTime[1]
    tEnd = srhhydro_SimTime[2] * 3600 #convert hours to seconds
    dt = srhhydro_SimTime[3]

    #update the constants
    swe_2D_constants.tStart = tStart
    swe_2D_constants.tEnd = tEnd
    swe_2D_constants.dt = dt
    swe_2D_constants.t = tStart   #set the current time as the starting time 
end