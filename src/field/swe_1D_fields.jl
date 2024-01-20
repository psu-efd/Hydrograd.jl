
#Variables and solution fields for 1D SWE with finite volume method


Base.@kwdef mutable struct swe_1D_parameters
    g::Float64 = 9.81   #gravity constant 

    #ManningN     #Manning's n

    t::Float64   #time 
    dt_min::Float64 = 0.001  #minimum time step size 
    dt::Float64   #time step size 
    CFL::Float64 = 0.4  #CFL number 
    tEnd::Float64  #end time for simulation 

    h_small::Float64 = 1.0e-3  #dry bed water depth threshold (e.g., 1.0e-3)

    RiemannSolver::String = "HLL"  #The choice of Riemann solver, e.g., Roe, HLL, and HLLC 
end

#enum of available boundary type names 
#@enum Boundary_Type_Name begin
#    internal   #not even an external boundary 
#    wall
#    zeroGradient
#    inletQ
#    exitH
#end

#struct for a boundary 
# bcType, bcValue:
    #    0-wall (no flux, dh/dx=0 and q=hu=0), bcValue = 0.0 (not used)
    #    1-zeroGradient (dh/dx=0, dq/dx=0), bcValue = 0.0 (not used)
    #    2-inletQ (specify inlet specific discharge q=hu=q_in), bcValue = q_in
    #    3-exitH (specify outlet water depth, h = h_outlet), bcValue = h_outlet
#Base.@kwdef mutable struct Boundary_1D
#    bcType::String
#    bcValue
#end

Base.@kwdef mutable struct swe_1D_fields{T<:Number}  
    
    h::Vector{T}    # water depth
    eta::Vector{T}  # free surface elevation: eta = h + zb
    q::Vector{T}    # the conservative variable hu

    zb_face::Vector{T}    # bed elevation at faces (points in 1D)
    zb_cell::Vector{T}    # bed elevation at cell centers
    S0::Vector{T}    # bed slope at cell centers 
    
    total_water_volume::Vector{T}  #total water volume in the domain 
end

# Function to initialize the swe_1D_fields struct
function initialize_swe_1D_fields(mesh::mesh_1D)
    
    nCells = mesh.nCells
    nFaces = mesh.nFaces

    T = typeof(1.0)

    return swe_1D_fields(
                    h=zeros(T,nCells), eta=zeros(T,nCells), q=zeros(T,nCells), 
                    zb_face=zeros(T,nFaces), zb_cell=zeros(T,nCells), S0=zeros(T,nCells),
                    total_water_volume = Float64[])
end

