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

    #input_dim = 1
    #output_dim = 1
    #hidden_layers = [3, 3]
    #activations = ["leakyrelu", "leakyrelu"]

    # Construct the model
    layers = []
    in_dim = input_dim
    for (i, hidden_dim) in enumerate(hidden_layers)
        push!(layers, Dense(in_dim, hidden_dim, get_activation(activations[i]); init_weight = Lux.glorot_uniform, init_bias = Lux.zeros64))
        push!(layers, LayerNorm(hidden_dim))  # Add normalization
        in_dim = hidden_dim
    end
    push!(layers, Dense(in_dim, output_dim; init_weight = Lux.glorot_uniform, init_bias = Lux.zeros64))  # Output layer without activation

     # Wrap model to enforce output bounds
     UDE_model = Chain(layers..., (x -> output_bounds[1] .+ (output_bounds[2] - output_bounds[1]) .* sigmoid(x)))

     # Randomly initialize the NN model
     UDE_model_params, UDE_model_state = Lux.setup(rng, UDE_model)

     print_lux_network_structure(UDE_model, UDE_model_params)

    if settings.UDE_settings.UDE_NN_config["how_to_initialize_NN"] == "random"
        # do nothing (it is already randomly initialized above)        
    elseif settings.UDE_settings.UDE_NN_config["how_to_initialize_NN"] == "from_pretrained"
        UDE_model_params_loaded, UDE_model_state_loaded = load_NN_model(settings.UDE_settings.UDE_NN_config["NN_weights_state_file_name"])

        #check the compatibility of the NN model with the input and output dimensions
        if typeof(UDE_model_params) !== typeof(UDE_model_params_loaded) || # Compare types
            length(UDE_model_params) !== length(UDE_model_params_loaded)  # Compare lengths
            error("The loaded pretrained NN model parameters and state do not match the model setup.")
        end

        #copy the loaded parameters and state to the model
        UDE_model_params = UDE_model_params_loaded
        UDE_model_state = UDE_model_state_loaded
     else
        error("Invalid how_to_initialize_NN: $(settings.UDE_settings.UDE_NN_config["how_to_initialize_NN"]). Supported options: random, from_pretrained.")
     end

     return UDE_model, UDE_model_params, UDE_model_state
end

function print_lux_network_structure(model, UDE_model_params)
    println("\nLux Neural Network Structure:")
    println("-------------------------")
    
    if model isa Lux.Chain
        for (i, layer) in enumerate(model.layers)
            println("Layer $i: ", typeof(layer))
            if layer isa Lux.Dense
                println("  Input size:  ", layer.in_dims)
                println("  Output size: ", layer.out_dims)
                println("  Activation:  ", layer.activation)
            end
        end
    end
    
    # Print parameter shapes
    println("\nParameter shapes:")
    for (k, v) in pairs(UDE_model_params)
        println("  $k: ", v)
    end
    println("-------------------------")
end

# Function to map activation function names to Lux activation functions
function get_activation(name::String)
    if name == "relu"
        return relu
    elseif name == "leakyrelu"
        return leakyrelu
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
function save_NN_model(model, params, state, filename::String)
    jldsave(filename; model, params, state)
    println("   NN model and parameters saved to $filename")
end

# Load model and parameters
function load_NN_model(filename::String)
    data = jldopen(filename)
    return data[:params], data[:state]
end

# Example usage (after training):
# Assuming `params` contains trained parameters
# save_model(model, params, "trained_model.bson")

# To reload:
# loaded_model, loaded_params = load_model("trained_model.bson")
# println("Loaded model: $loaded_model")
