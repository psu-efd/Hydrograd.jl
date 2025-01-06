#process the Manning's n values 

function setup_ManningN(settings, my_mesh_2D, srh_all_Dict)

    # Initialize Manning's n values for cells
    ManningN_cells = zeros(Float64, my_mesh_2D.numOfCells)

    srhhydro_ManningsN = srh_all_Dict["srhhydro_ManningsN"]
    srhmat_numOfMaterials = srh_all_Dict["srhmat_numOfMaterials"]
    srhmat_matNameList = srh_all_Dict["srhmat_matNameList"]
    srhmat_matZoneCells = srh_all_Dict["srhmat_matZoneCells"]
    matID_cells = srh_all_Dict["matID_cells"]

    if settings.bVerbose
        println("srhhydro_ManningsN: ", srhhydro_ManningsN)
        println("srhmat_numOfMaterials: ", srhmat_numOfMaterials)
        println("srhmat_matNameList: ", srhmat_matNameList)
        println("srhmat_matZoneCells: ", srhmat_matZoneCells)
        println("matID_cells: ", matID_cells)
    end

    #loop over all cells to setup Manning's n
    for iCell in 1:my_mesh_2D.numOfCells
        matID = matID_cells[iCell]
        matName = srhmat_matNameList[matID]

        if haskey(srhhydro_ManningsN, matID)
            ManningN_cells[iCell] = srhhydro_ManningsN[matID]
        else
            error("Material $matName does not have Manning's n. Assign the default value 0.03.")
        end
    end

    #println("ManningN_cells: ", ManningN_cells)

    return ManningN_cells
end

#update Manning's n values for inversion or sensitivity analysis
# UPdate based on the provided Manning's n values for each material (zone)
# new_ManningN_values is a vector of Manning's n values for each material (zone)
function update_ManningN_inversion_sensitivity_analysis(
    my_mesh_2D::mesh_2D,
    srh_all_Dict::Dict{String,Any},
    new_ManningN_values::Vector{T}) where {T<:Real}


    #material ID for each cell (0-based): 0-default material, 1-first material, 2-second material, etc.
    matID_cells = srh_all_Dict["matID_cells"]

    # Create array directly with comprehension
    ManningN_cells = [new_ManningN_values[matID_cells[i]+1] for i in 1:my_mesh_2D.numOfCells]  #+1 to make matID_cells 1-based (to be consistent with the new_ManningN_values)


    return ManningN_cells
end

#update Manning's n values based on the provided Manning's n values for each material (zone)
# new_ManningN_values is a vector of Manning's n values for each material (zone)
function update_ManningN_forward_simulation(
    h::AbstractArray{T},
    settings::ControlSettings) where {T<:Real}

    manning_func = create_manning_function(
            settings.forward_settings.ManningN_function_type,
            settings.forward_settings.ManningN_function_parameters
        )

    ManningN_cells = manning_func(h)

    return ManningN_cells
end

# Function factory to create Manning's n function from parameters
function create_manning_function(type::String, params::Dict)
    if type == "power_law"
        return h -> power_law_n(h, Float64(params["n_lower"]), Float64(params["n_upper"]), Float64(params["k"]))
    elseif type == "sigmoid"
        return h -> sigmoid_n(h, Float64(params["n_lower"]), Float64(params["n_upper"]), Float64(params["k"]), Float64(params["h_mid"]))
    elseif type == "inverse"
        return h -> inverse_n(h, Float64(params["n_lower"]), Float64(params["n_upper"]), Float64(params["k"]))
    else
        error("Unknown Manning's n function type: $type. Supported types: constant, power_law, sigmoid, inverse.")
    end
end

# Predefined functions for Manning's n(h)
# Inverse function: n(h) = n_lower + (n_upper - n_lower) / (1 + k*h)
#    k is a positive scaling factor that contrrols the the rate of decrease of Manning's n with increasing water depth
function inverse_n(h::AbstractArray{T}, n_lower::T, n_upper::T, k::T) where {T<:Real}
    @assert(k > 0, "k must be a positive scaling factor")
    return n_lower .+ (n_upper - n_lower) ./ (1 .+ k .* h)
end

# Power law function: n(h) = n_lower + (n_upper - n_lower) * h^(-k)
#    k is a positive scaling factor that contrrols the the rate of decrease of Manning's n with increasing water depth
function power_law_n(h::AbstractArray{T}, n_lower::T, n_upper::T, k::T) where {T<:Real}
    @assert(k > 0, "k must be a positive scaling factor")
    return n_lower .+ (n_upper - n_lower) .* (h .+ eps(eltype(h))) .^ (-k)
end

# sigmoid based on water depth: n(h) = n_lower + (n_upper - n_lower) / (1 + exp(-k*(h-h_mid)))
#    k is a positive scaling factor that contrrols the the rate of decrease of Manning's n with increasing water depth
#    h_mid is the water depth at which Manning's n transitions from n_lower to n_upper rapidly.
function sigmoid_n(h::AbstractArray{T}, n_lower::T, n_upper::T, k::T, h_mid::T) where {T<:Real}
    @assert(k > 0, "k must be a positive scaling factor")
    @assert(h_mid > 0, "h_mid must be a positive water depth")
    return n_lower .+ (n_upper - n_lower) ./ (1 .+ exp.(k .* (h .- h_mid)))
end

#update Manning's n values from the UDE model's neural network
function update_ManningN_UDE(h, ude_model, NN_model_params, ude_model_state, num_of_cells)
    # Need to convert scalar input (h) to matrices for Lux
    ManningN_cells = [
        let
            # Reshape scalar h[iCell] into a 1×1 matrix
            h_matrix = reshape([h[iCell]], 1, 1)
            # Apply model and extract scalar result
            ude_model(h_matrix, NN_model_params, ude_model_state)[1][1]
        end
        for iCell in 1:num_of_cells
    ]

    Zygote.ignore() do
        #@show typeof(ManningN_cells)
        #@show size(ManningN_cells)
        #@show ManningN_cells
    end

    return ManningN_cells
end