
#update Flow resistance values from the UDE model's neural network: Flow resistance = f(h, u, v)
function update_FlowResistance_UDE(h::AbstractVector{T},
    q_x::AbstractVector{T},
    q_y::AbstractVector{T},
    ude_model::Lux.Chain,
    NN_model_params,
    ude_model_state,
    h_bounds::Vector{T},
    velocity_magnitude_bounds::Vector{T},
    num_of_cells::Integer) where {T<:Real}

    #compute the bounds for q
    q_bounds = [h_bounds[1]*velocity_magnitude_bounds[1], h_bounds[2]*velocity_magnitude_bounds[2]]

    #compute the magnitude of q
    q_magnitude = smooth_sqrt.(q_x.^2 + q_y.^2)

    # Create batch of inputs and normalize to [-1, 1]
    h_normalized = @. 2.0 * (h - h_bounds[1]) / (h_bounds[2] - h_bounds[1]) - 1.0
    q_magnitude_normalized = @. 2.0 * (q_magnitude - q_bounds[1]) / (q_bounds[2] - q_bounds[1]) - 1.0

    # Reshape for Lux
    h_matrix = reshape(h_normalized, 1, num_of_cells)
    q_magnitude_matrix = reshape(q_magnitude_normalized, 1, num_of_cells)
    
    #concatenate the input variables
    input_matrix = vcat(h_matrix, q_magnitude_matrix)

    # Apply model to entire batch
    outputs = ude_model(input_matrix, NN_model_params, ude_model_state)[1]
    
    friction_magnitudes = vec(outputs)

    Zygote.ignore() do
        #@show typeof(input_matrix)
        #@show size(input_matrix)
        #@show input_matrix

        #@show typeof(friction_magnitudes)
        #@show size(friction_magnitudes)
        #@show friction_magnitudes
    end
    
    return friction_magnitudes
end