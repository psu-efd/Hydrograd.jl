using JLD2
using Lux
using Zygote
using StableRNGs
using JSON3
using ComponentArrays

using LaTeXStrings  
using Plots

gr()

#default plot setup 
plot_font = "Computer Modern"
#plot_font = "times"
default(
    fontfamily=plot_font,
    linewidth=2, 
    framestyle=:box, 
    label=nothing, 
    grid=false
)

# Scale all font sizes by a factor
Plots.scalefontsizes()
Plots.scalefontsizes(10.0)


#compare the resulted n(h) from the trained NN with the truth


#compute the truth n(h)
function truth_n_h(h_values)

    #truth n(h) is a sigmoid function of h (this is hardcoded here)
    truth_n_h_parameters = Dict(
        "n_lower" => 0.03,
        "n_upper" => 0.06,
        "k" => 100.0,
        "h_mid" => 0.3
    )

    #compute the truth n(h)
    truth_n_h = truth_n_h_parameters["n_lower"] .+ (truth_n_h_parameters["n_upper"] - truth_n_h_parameters["n_lower"]) ./ (1 .+ exp.(truth_n_h_parameters["k"] .* (h_values .- truth_n_h_parameters["h_mid"])))

    return truth_n_h
end

#trained n(h) is from the NN
function NN_n_h(h_values)
    # read the control file
    control_file = "run_control.json"
    control_dict = JSON3.read(open(control_file), Dict)

    #load the UDE training results 
    #UDE_training_results_file_name = control_dict["UDE_options"]["UDE_NN_config"]["NN_weights_state_file_name"]
    UDE_training_results_file_name = "UDE_training_results.jld2"
    UDE_training_results = load(UDE_training_results_file_name)

    #Load the UDE model, parameters (only the last iteration), and state
    UDE_model_params_pretrained = UDE_training_results["ude_model_params"]
    UDE_model_state_pretrained = UDE_training_results["ude_model_state"]

    # Set a random seed for reproducible behaviour
    rng = StableRNG(1111)

    # NN configuration
    input_dim = control_dict["UDE_options"]["UDE_NN_config"]["input_dim"]
    output_dim = control_dict["UDE_options"]["UDE_NN_config"]["output_dim"]
    hidden_layers = control_dict["UDE_options"]["UDE_NN_config"]["hidden_layers"]
    activations = control_dict["UDE_options"]["UDE_NN_config"]["activations"]
    output_bounds = control_dict["UDE_options"]["UDE_NN_config"]["output_bounds"]

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

    # Randomly initialize the NN model
    UDE_model_params, UDE_model_state = Lux.setup(rng, UDE_model)

    #use the pretrained parameters and state to make the NN prediction
    n_h_NN = update_ManningN_UDE(h_values, UDE_model, UDE_model_params_pretrained, UDE_model_state_pretrained, length(h_values))

    return n_h_NN
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

#plot and compare the truth and the trained n(h)
function plot_and_compare_n_h(h_values, truth_n_h_values, NN_n_h_values)

    # Create the plot
    p = plot(truth_n_h_values, h_values, label="Truth", linewidth=2, linestyle=:solid, color=:black,
			size=(800, 600),  # Figure size: 800x600 pixels
			dpi=300,          # High-resolution
			background_color=:white,  # Background color
			tickfontsize=12,    # Customize tick font size
			guidefontsize=14    # Customize labels' font size
	)
    plot!(p, NN_n_h_values, h_values, label="NN", linewidth=2, linestyle=:dash, color=:red)

    # Add labels and title
    xlabel!(p, L"Manning's $n$", fontsize=14)
    ylabel!(p, L"$h$ (m)", fontsize=14)
    #title!(p, "Comparsion of Truth and NN Predictions", fontsize=16)

	# Customize legend
    plot!(
        p, legend=:topright,  # Position legend at the top-right
		legendfontsize=12,  # Customize legend font size
        #legendtitle="Legend", # Add a legend title
        #legendfont=:sans      # Use a sans-serif font for the legend
        foreground_color_legend = nothing
    )
	
	# Modify the axis ranges
	xlims!(p, 0.0, 0.061) 
	ylims!(p, 0.0, 1.0)   

    # Save the plot
    savefig("compare_ManningN_h_Truth_NN.png")  # Save the plot as a PNG file

    # Display the plot
    display(p)
end


#define the h values
h_values = collect(range(0.001, 1.0, length=100))

#get the truth n(h)
truth_n_h_values = truth_n_h(h_values)

#get the trained n(h)
NN_n_h_values = NN_n_h(h_values)

@show size(h_values)
@show h_values
@show size(truth_n_h_values)
@show truth_n_h_values
@show size(NN_n_h_values)
@show NN_n_h_values

#plot the results
plot_and_compare_n_h(h_values, truth_n_h_values, NN_n_h_values)
