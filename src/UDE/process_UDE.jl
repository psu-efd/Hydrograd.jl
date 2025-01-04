using Lux, Random, JSON3

# Function to create a neural network from configuration: UDE_options["UDE_NN_config"]
function create_NN_model(settings)

    if !settings.bPerform_UDE
        return nothing, nothing, nothing
    end

    # Set a random seed for reproducible behaviour
    rng = StableRNG(1111)

    # Extract configuration settings
    input_dim = settings.UDE_settings.UDE_NN_config["input_dim"]
    output_dim = settings.UDE_settings.UDE_NN_config["output_dim"]
    hidden_layers = settings.UDE_settings.UDE_NN_config["hidden_layers"]  # List of hidden layer sizes
    activations = settings.UDE_settings.UDE_NN_config["activations"]  # List of activation functions
    output_bounds = settings.UDE_settings.UDE_NN_config["output_bounds"]  # List of output bounds

    # Validate configuration
    if length(hidden_layers) != length(activations)
        error("The number of hidden layers must match the number of activation functions.")
    end

    # Construct the model
    layers = []
    in_dim = input_dim
    for (i, h_dim) in enumerate(hidden_layers)
        push!(layers, Dense(in_dim, h_dim, get_activation(activations[i])))
        in_dim = h_dim
    end
    push!(layers, Dense(in_dim, output_dim))  # Output layer without activation

     # Wrap model to enforce output bounds
     UDE_model = Chain(layers..., (x -> output_bounds[1] .+ (output_bounds[2] - output_bounds[1]) .* sigmoid(x)))

     UDE_model_params, UDE_model_state = Lux.setup(rng, UDE_model)

     return UDE_model, UDE_model_params, UDE_model_state
end

# Function to map activation function names to Lux activation functions
function get_activation(name::String)
    if name == "relu"
        return relu
    elseif name == "sigmoid"
        return sigmoid
    elseif name == "tanh"
        return tanh
    elseif name == "softplus"
        return softplus
    else
        error("Unsupported activation function: $name")
    end
end

# Save model and parameters
function save_NN_model(model, params, filename::String)
#    BSON.@save filename model params
    println("Model and parameters saved to $filename")
end

# Load model and parameters
function load_NN_model(filename::String)
#    data = BSON.@load filename
    return data[:model], data[:params]
end

# Example usage (after training):
# Assuming `params` contains trained parameters
# save_model(model, params, "trained_model.bson")

# To reload:
# loaded_model, loaded_params = load_model("trained_model.bson")
# println("Loaded model: $loaded_model")
