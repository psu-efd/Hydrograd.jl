#process the Manning's n values 

function setup_ManningN(settings, my_mesh_2D, srh_all_Dict)

    # Initialize Manning's n values for cells (just for initialization; will be updated later if Manning's n is a function of h, ks, or Umag)
    ManningN_cells = zeros(Float64, my_mesh_2D.numOfCells)
    ks_cells = zeros(Float64, my_mesh_2D.numOfCells)
    h_ks_cells = zeros(Float64, my_mesh_2D.numOfCells)
    friction_factor_cells = zeros(Float64, my_mesh_2D.numOfCells)
    Re_cells = zeros(Float64, my_mesh_2D.numOfCells)


    srhhydro_ManningsN = srh_all_Dict["srhhydro_ManningsN"]
    srhmat_numOfMaterials = srh_all_Dict["srhmat_numOfMaterials"]
    srhmat_matNameList = srh_all_Dict["srhmat_matNameList"]
    srhmat_matZoneCells = srh_all_Dict["srhmat_matZoneCells"]
    srhmat_ks_dict = srh_all_Dict["srhmat_ks_dict"]    #this is extra roughness height information for each material (not in srhmat file)
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

    #if it is forward simulation and the Manning's n variable, or for UDE and UDE_choice is ManningN_h_Umag_ks, then create the variable ks_cells
    bNeed_ks_cells = false
    if settings.bPerform_Forward_Simulation 
        if settings.forward_settings.ManningN_option == "variable" 
            bNeed_ks_cells = true
        end
    end

    if settings.bPerform_UDE 
        if settings.UDE_settings.UDE_choice == "ManningN_h_Umag_ks" 
            bNeed_ks_cells = true
        end
    end

    println("srhmat_ks_dict: ", srhmat_ks_dict)

    if bNeed_ks_cells
        for iCell in 1:my_mesh_2D.numOfCells
            matID = matID_cells[iCell]
            ks_cells[iCell] = srhmat_ks_dict[matID]
        end        
    end

    #println("ManningN_cells: ", ManningN_cells)

    return ManningN_cells, ks_cells, h_ks_cells, friction_factor_cells, Re_cells
end

#update Manning's n values for inversion or sensitivity analysis
# Update based on the provided Manning's n values for each material (zone)
# new_ManningN_values is a vector of Manning's n values for each material (zone)
function update_ManningN_inversion_sensitivity_analysis(
    my_mesh_2D::mesh_2D,
    srh_all_Dict::Dict{String,Any},
    new_ManningN_values::Vector{T}) where {T<:Real}


    #material ID for each cell (0-based): 0-default material, 1-first material, 2-second material, etc.
    matID_cells = srh_all_Dict["matID_cells"]

    Zygote.ignore() do
        #println("In update_ManningN_inversion_sensitivity_analysis:")
        #println("  new_ManningN_values size: ", size(new_ManningN_values))
        #println("  num cells: ", my_mesh_2D.numOfCells)
        #println("  matID_cells size: ", size(srh_all_Dict["matID_cells"]))
    end

    # Create array directly with comprehension
    ManningN_cells = [new_ManningN_values[matID_cells[i]+1] for i in 1:my_mesh_2D.numOfCells]  #+1 to make matID_cells 1-based (to be consistent with the new_ManningN_values)

    Zygote.ignore() do
        #println("  ManningN_cells size: ", size(ManningN_cells))
        #println("  Sample values:")
        #println("    First 3 matID_cells: ", matID_cells[1:3])
        #println("    First 3 ManningN_cells: ", ManningN_cells[1:3])
    end

    return ManningN_cells
end

#update Manning's n values based on the provided Manning's n values for each material (zone)
# new_ManningN_values is a vector of Manning's n values for each material (zone)
function update_ManningN_forward_simulation(
    h::AbstractArray{T},
    Umag::AbstractArray{T},
    ks::AbstractArray{T},
    settings::ControlSettings) where {T<:Real}

    manning_func = create_manning_function(
            settings.forward_settings.ManningN_function_type,
            settings.forward_settings.ManningN_function_parameters
        )

    #ManningN_cells = manning_func(h, Umag, ks)
    ManningN_cells, h_ks, friction_factor, Re = manning_func(h, Umag, ks)   #debug

    return ManningN_cells, h_ks, friction_factor, Re
end

# Function factory to create Manning's n function from parameters
# For types "power_law", "sigmoid", and "inverse", the function is a function of h only (Umag and ks are not used)
# For type "h_Umag_ks", the function is a function of h, Umag, and ks
function create_manning_function(type::String, params::Dict)
    if type == "power_law"
        return (h, Umag, ks) -> power_law_n(h, Float64(params["n_lower"]), Float64(params["n_upper"]), Float64(params["k"]))
    elseif type == "sigmoid"
        return (h, Umag, ks) -> sigmoid_n(h, Float64(params["n_lower"]), Float64(params["n_upper"]), Float64(params["k"]), Float64(params["h_mid"]))
    elseif type == "inverse"
        return (h, Umag, ks) -> inverse_n(h, Float64(params["n_lower"]), Float64(params["n_upper"]), Float64(params["k"]))
    elseif type == "h_Umag_ks"
        return (h, Umag, ks) -> manning_n_h_Umag_ks(h, Umag, ks)
    else
        error("Unknown Manning's n function type: $type. Supported types: constant, power_law, sigmoid, inverse.")
    end
end

# Predefined functions for Manning's n(h)
# Inverse function: n(h) = n_lower + (n_upper - n_lower) / (1 + k*h)
#    k is a positive scaling factor that contrrols the the rate of decrease of Manning's n with increasing water depth
function inverse_n(h::AbstractArray{T}, n_lower::T, n_upper::T, k::T) where {T<:Real}
    @assert(k > 0, "k must be a positive scaling factor")
    
    #h_ks, f, and Re are not used for inverse_n
    h_ks = zeros(eltype(h), size(h))
    f = zeros(eltype(h), size(h))
    Re = zeros(eltype(h), size(h))

    n = n_lower .+ (n_upper - n_lower) ./ (1 .+ k .* h)

    return n, h_ks, f, Re
end

# Power law function: n(h) = n_lower + (n_upper - n_lower) * h^(-k)
#    k is a positive scaling factor that contrrols the the rate of decrease of Manning's n with increasing water depth
function power_law_n(h::AbstractArray{T}, n_lower::T, n_upper::T, k::T) where {T<:Real}
    @assert(k > 0, "k must be a positive scaling factor")

    #h_ks, f, and Re are not used for power_law_n
    h_ks = zeros(eltype(h), size(h))
    f = zeros(eltype(h), size(h))
    Re = zeros(eltype(h), size(h))

    n = n_lower .+ (n_upper - n_lower) .* (h .+ eps(eltype(h))) .^ (-k)

    return n, h_ks, f, Re
end

# sigmoid based on water depth: n(h) = n_lower + (n_upper - n_lower) / (1 + exp(-k*(h-h_mid)))
#    k is a positive scaling factor that contrrols the the rate of decrease of Manning's n with increasing water depth
#    h_mid is the water depth at which Manning's n transitions from n_lower to n_upper rapidly.
function sigmoid_n(h::AbstractArray{T}, n_lower::T, n_upper::T, k::T, h_mid::T) where {T<:Real}
    @assert(k > 0, "k must be a positive scaling factor")
    @assert(h_mid > 0, "h_mid must be a positive water depth")

    #h_ks, f, and Re are not used for sigmoid_n
    h_ks = zeros(eltype(h), size(h))
    f = zeros(eltype(h), size(h))
    Re = zeros(eltype(h), size(h))

    n = n_lower .+ (n_upper - n_lower) ./ (1 .+ exp.(k .* (h .- h_mid)))

    return n, h_ks, f, Re
end

# Manning's n as a function of h, ks, and Umag
# See Cheng (2008) JHE, "Formulas for Friction Factor in Transitional Regimes", Eq. (17) to get f, then the Manning's n
function manning_n_h_Umag_ks(h::AbstractArray{T}, Umag::AbstractArray{T}, ks::AbstractArray{T}) where {T<:Real}

    # Computer the Reynolds number
    ν = 1.0e-6 # kinematic viscosity of water (m^2/s)
    Re = Umag .* h ./ ν

    # Compute h/ks 
    h_ks = h ./ ks

    #compute alpha
    alpha = 1.0 ./ (1.0 .+ (Re./850.0).^9)

    #compute beta
    beta = 1.0 ./ (1.0 .+ (Re./(h_ks.*160.0)).^2)

    # Compute the friction factor f in parts 
    part1 = (Re./24.0).^alpha 
    part2 = (1.8 .* log10.(Re./2.1)).^(2.0.*(1.0 .- alpha).*beta)
    part3 = (2.0 .* log10.(11.8 .* h_ks)).^(2.0.*(1.0 .- alpha).*(1.0 .- beta))

    # Compute the friction factor f
    f = 1.0 ./ (part1 .* part2 .* part3)

    # Compute Manning's n = sqrt(f/8.0) * h^(1/6) /sqrt(9.81)
    n = sqrt.(f./8.0) .* h.^(1.0./6.0) ./ sqrt.(9.81)       #need to fix this if in English units

    return n, h_ks, f, Re
end

#update Manning's n values from the UDE model's neural network
function update_ManningN_UDE(UDE_choice::String,
                           h::AbstractVector{T1}, 
                           Umag::AbstractVector{T2},
                           ks::AbstractVector{T3},
                           ude_model::Lux.Chain, 
                           NN_model_params::AbstractVector{T4}, 
                           ude_model_state,
                           h_bounds::Vector{Float64},
                           Umag_bounds::Vector{Float64},
                           ks_bounds::Vector{Float64},
                           num_of_cells::Integer) where {T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real}    

    #ManningN_cells = [
    #    ude_model(@SMatrix([h[iCell];;]), NN_model_params, ude_model_state)[1][1]
    #    for iCell in 1:num_of_cells
    #]

    # Create batch of inputs and normalize to [-1, 1]
    h_normalized = @. 2.0 * (h - h_bounds[1]) / (h_bounds[2] - h_bounds[1]) - 1.0

    if UDE_choice == "ManningN_h_Umag_ks"
        Umag_normalized = @. 2.0 * (Umag - Umag_bounds[1]) / (Umag_bounds[2] - Umag_bounds[1]) - 1.0
        ks_normalized = @. 2.0 * (ks - ks_bounds[1]) / (ks_bounds[2] - ks_bounds[1]) - 1.0
    end

    # Reshape for Lux
    h_matrix = reshape(h_normalized, 1, num_of_cells)

    if UDE_choice == "ManningN_h_Umag_ks"
        Umag_matrix = reshape(Umag_normalized, 1, num_of_cells)
        ks_matrix = reshape(ks_normalized, 1, num_of_cells)
    end

    Zygote.ignore() do
        #print a few values
        #println("h_matrix: ", ForwardDiff.value.(h_matrix[1:5]))
        #println("Umag_normalized: ", ForwardDiff.value.(Umag_normalized[1:5]))
        #println("ks_normalized: ", ForwardDiff.value.(ks_normalized[1:5]))
    end

    # Apply model to entire batch
    if UDE_choice == "ManningN_h"
        outputs = ude_model(h_matrix, NN_model_params, ude_model_state)[1]
    elseif UDE_choice == "ManningN_h_Umag_ks"
        outputs = ude_model(vcat(h_normalized', Umag_normalized', ks_normalized'), NN_model_params, ude_model_state)[1]
    else
        error("Unknown UDE choice: $UDE_choice")
    end

    Zygote.ignore() do
        #@show typeof(ManningN_cells)
        #@show size(ManningN_cells)
        #@show ManningN_cells
    end

    #return ManningN_cells
    return vec(outputs)
end