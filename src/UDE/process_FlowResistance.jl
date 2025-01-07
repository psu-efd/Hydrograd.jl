
#update Flow resistance values from the UDE model's neural network: Flow resistance = f(h, u, v)
function update_FlowResistance_UDE(h::AbstractVector{T}, 
    u::AbstractVector{T},
    v::AbstractVector{T},
    ude_model::Lux.Chain, 
    NN_model_params, 
    ude_model_state,
    h_bounds::Vector{T},
    u_bounds::Vector{T},
    v_bounds::Vector{T},
    num_of_cells::Integer) where T <: Real


# Create batch of inputs and normalize to [-1, 1]
h_normalized = @. 2.0 * (h - h_bounds[1]) / (h_bounds[2] - h_bounds[1]) - 1.0
u_normalized = @. 2.0 * (u - u_bounds[1]) / (u_bounds[2] - u_bounds[1]) - 1.0
v_normalized = @. 2.0 * (v - v_bounds[1]) / (v_bounds[2] - v_bounds[1]) - 1.0

# Reshape for Lux
h_matrix = reshape(h_normalized, 1, num_of_cells)
u_matrix = reshape(u_normalized, 1, num_of_cells)
v_matrix = reshape(v_normalized, 1, num_of_cells)

input_matrix = vcat(h_matrix, u_matrix, v_matrix)

# Apply model to entire batch
outputs = ude_model(input_matrix, NN_model_params, ude_model_state)[1]

Zygote.ignore() do
#@show typeof(ManningN_cells)
#@show size(ManningN_cells)
#@show ManningN_cells
end

#return ManningN_cells
return vec(outputs)
end