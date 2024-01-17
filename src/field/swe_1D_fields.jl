
#Variables and solution fields for 1D SWE with finite volume method


Base.@kwdef mutable struct swe_1D_parameters
    g::Float64 = 9.81   #gravity constant 

    ManningN::Float64 = 0.03    #Manning's n

    t::Float64   #time 
    dt_min::Float64 = 0.001  #minimum time step size 
    dt::Float64   #time step size 
    CFL::Float64 = 0.4  #CFL number 
    tEnd::Float64  #end time for simulation 

    h_small::Float64 = 1.0e-3  #dry bed water depth threshold (e.g., 1.0e-3)

    RiemannSolver::String = "HLL"  #The choice of Riemann solver, e.g., Roe, HLL, and HLLC 
end

#enum of available boundary type names 
@enum Boundary_Type_Name begin
    internal   #not even an external boundary 
    wall
    zeroGradient
    inletQ
    exitH
end

#struct for a boundary 
# bcType, bcValue:
    #    0-wall (no flux, dh/dx=0 and q=hu=0), bcValue = 0.0 (not used)
    #    1-zeroGradient (dh/dx=0, dq/dx=0), bcValue = 0.0 (not used)
    #    2-inletQ (specify inlet specific discharge q=hu=q_in), bcValue = q_in
    #    3-exitH (specify outlet water depth, h = h_outlet), bcValue = h_outlet
Base.@kwdef mutable struct Boundary_1D 
    bcType::Boundary_Type_Name
    bcValue::Float64
end

Base.@kwdef mutable struct swe_1D_fields  
    
    h::Vector{Float64}    # water depth
    eta::Vector{Float64}  # free surface elevation: eta = h + zb
    q::Vector{Float64}    # the conservative variable hu

    zb_face::Vector{Float64}    # bed elevation at faces (points in 1D)
    zb_cell::Vector{Float64}    # bed elevation at cell centers
    S0::Vector{Float64}    # bed slope at cell centers 
    
    dhdt::Vector{Float64}  #dh/dt 
    dqdt::Vector{Float64}  #dq/dt 

    h_old::Vector{Float64}  #h at previous time step/iteration 
    q_old::Vector{Float64}  #q at previous time step/iteration 
    dhdt_old::Vector{Float64}  #dh/dt at previous time step/iteration 
    dqdt_old::Vector{Float64}  #dq/dt at previous time step/iteration 

    fluxes::Array{Float64}  #flux, a 2D array with the shape of [2, nFaces]. flux[1, iFace] and flux[2, iFace] are for fluxes of cty and mom-x.

    total_water_volume::Vector{Float64}  #total water volume in the domain 

    #parameters
    swe_1d_para::swe_1D_parameters

    #boundary conditions: for 1D, we only have left and right boundaries 
    leftBoundary::Boundary_1D
    rightBoundary::Boundary_1D    
end

# Function to initialize the swe_1D_fields struct
function initialize_swe_1D_fields(mesh::mesh_1D, swe_1d_para::swe_1D_parameters, 
    leftBoundary::Boundary_1D, rightBoundary::Boundary_1D)
    
    nCells = mesh.nCells
    nFaces = mesh.nFaces

    return swe_1D_fields(
                    h=zeros(nCells), eta=zeros(nCells), q=zeros(nCells), zb_face=zeros(nFaces), zb_cell=zeros(nCells), S0=zeros(nCells),
                    dhdt=zeros(nCells), dqdt=zeros(nCells), h_old=zeros(nCells), q_old=zeros(nCells), dhdt_old=zeros(nCells), dqdt_old=zeros(nCells),
                    fluxes=zeros(2,nFaces),
                    total_water_volume = Float64[], swe_1d_para=swe_1d_para, leftBoundary=leftBoundary, rightBoundary=rightBoundary )
end

