mutable struct BoundaryConditions
    nInletQ::Int
    inletQ_BC_indices::Vector{Int}
    inletQ_faceCentroids::Vector{Matrix{Float64}}

    ## Inner constructor with keyword arguments
    function BoundaryConditions(;
        nInletQ::Int=0,
        inletQ_BC_indices::Vector{Int}=Int[],
        inletQ_faceCentroids::Vector{Matrix{Float64}}=Matrix{Float64}[]
    )
        return new(nInletQ, inletQ_BC_indices, inletQ_faceCentroids)
    end

end

function create_boundary_conditions(srh_all_Dict)
    # Example fields
    nInletQ = srh_all_Dict["nInletQ_BCs"]
    inletQ_BC_indices = Int[]  # Initialize empty vector
    inletQ_faceCentroids = Vector{Matrix{Float64}}(undef, nInletQ)  # Pre-allocate

    #create an empty BoundaryConditions object
    boundary_conditions = BoundaryConditions()

    #modify the BoundaryConditions object
    boundary_conditions.nInletQ = nInletQ
    boundary_conditions.inletQ_BC_indices = inletQ_BC_indices
    boundary_conditions.inletQ_faceCentroids = inletQ_faceCentroids

    return boundary_conditions
end

srh_all_Dict = Dict("nInletQ_BCs" => 3)
boundary_conditions = create_boundary_conditions(srh_all_Dict)

println(boundary_conditions)
