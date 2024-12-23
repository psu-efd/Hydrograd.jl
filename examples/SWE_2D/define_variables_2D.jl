#define variables 

#mesh related variables
#nodeCoordinates::Array{Float64,2}  # Node coordinates: Float64 2D array [numOfNodes, 3]
nodeCoordinates = srhgeom_obj.nodeCoordinates

#cellBedElevation::Vector{Float64}  # Element bed elevation: Float64 1D array [numOfCells]
cellBedElevation = srhgeom_obj.elementBedElevation

#inlet-q boundary related variables
inletQ_TotalQ::Vector{Float64}                         #total discharge of the inlet-q boundaries
inletQ_H::Vector{Vector{Float64}}                     #water depth of the inlet-q boundaries
inletQ_A::Vector{Vector{Float64}}                     #cross-sectional area of the inlet-q boundaries
inletQ_ManningN::Vector{Vector{Float64}}             #Manning's n of the inlet-q boundaries
inletQ_Length::Vector{Vector{Float64}}               #length of the inlet-q boundaries
inletQ_TotalA::Vector{Float64}                       #total cross-sectional area of the inlet-q boundaries
inletQ_DryWet::Vector{Vector{Int}}                   #dry/wet flag of the inlet-q boundaries

#exit-h boundary related variables
exitH_WSE::Vector{Float64}                             #water surface elevation of the exit-h boundaries
exitH_H::Vector{Vector{Float64}}                     #water depth of the exit-h boundaries
exitH_A::Vector{Vector{Float64}}                     #cross-sectional area of the exit-h boundaries

#wall boundary related variables
wall_H::Vector{Vector{Float64}}                     #water depth of the wall boundaries
wall_A::Vector{Vector{Float64}}                     #cross-sectional area of the wall boundaries

#symmetry boundary related variables
symm_H::Vector{Vector{Float64}}                     #water depth of the symmetry boundaries
symm_A::Vector{Vector{Float64}}                     #cross-sectional area of the symmetry boundaries
