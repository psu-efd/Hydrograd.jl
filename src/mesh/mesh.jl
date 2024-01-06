#Mesh functions

using Base

#2D Mesh in SRHGEOM format
Base.@kwdef mutable struct Mesh2D_SRHGEOM
    Type::String = "SRHGEOM"
    Name::String = " "      #name of the mesh
    GridUnit::String = "SI" #unit of the mesh, SI or EN

    bcDict::Dict{Int64,String} = Dict{Int64,String}()             #boundary condition dictionary

    numOfElements::Integer = -1 #number of elements 
    numOfNodes::Integer = -1 #number of nodes  
    numOfNodeStrings::Integer = -1 #number of Nodestrings 

    elementNodesList::Vector{Int64} = [] #list of nodes for all elements 
end