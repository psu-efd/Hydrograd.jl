using Revise

using Hydrograd

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


#compare the resulted FlowResistance(h, u, v) from the trained NN with the truth

#compute the truth flow resistance
function truth_flow_resistance(h_values, u_values, v_values)
    # Example: Using Manning's formula
    # You can modify this with your actual truth function
    n = 0.03  # Manning's n
    g = 9.81  # gravity
    
    # Compute velocity magnitude
    u_mag = @. sqrt(u_values^2 + v_values^2 + eps())
    
    # Manning's formula for flow resistance
    tau_b = @. g * n^2 / (h_values^(1/3)) * u_mag

    return tau_b
end

#trained flow resistance from the NN
function NN_flow_resistance(h_values, u_values, v_values, settings)
    # Load the trained model
    UDE_training_results_file_name = "UDE_training_results.jld2"
    UDE_training_results = load(UDE_training_results_file_name)
    
    UDE_model_params_pretrained = UDE_training_results["ude_model_params"]
    UDE_model_state_pretrained = UDE_training_results["ude_model_state"]
    
    UDE_model, _, _ = Hydrograd.create_NN_model(settings)
    
    # Normalize inputs using bounds from settings
    h_bounds = Float64.(settings.UDE_settings.UDE_NN_config["h_bounds"])
    u_bounds = Float64.(settings.UDE_settings.UDE_NN_config["u_bounds"])
    v_bounds = Float64.(settings.UDE_settings.UDE_NN_config["v_bounds"])
    
    h_normalized = @. 2.0 * (h_values - h_bounds[1]) / (h_bounds[2] - h_bounds[1]) - 1.0
    u_normalized = @. 2.0 * (u_values - u_bounds[1]) / (u_bounds[2] - u_bounds[1]) - 1.0
    v_normalized = @. 2.0 * (v_values - v_bounds[1]) / (v_bounds[2] - v_bounds[1]) - 1.0
    
    # Reshape for NN input
    num_points = length(h_values)
    input_matrix = vcat(
        reshape(h_normalized, 1, num_points),
        reshape(u_normalized, 1, num_points),
        reshape(v_normalized, 1, num_points)
    )
    
    # Get NN prediction
    flow_resistance = vec(UDE_model(input_matrix, UDE_model_params_pretrained, UDE_model_state_pretrained)[1])
    
    return flow_resistance
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


#main 
case_path = pwd()

# Read and parse control file
println("Reading control file...")
control_file = joinpath(case_path, "run_control.json")
settings = Hydrograd.parse_control_file(control_file)


#define the h values
h_values = collect(range(0.001, 1.0, length=100))

#get the truth n(h)
truth_n_h_values = truth_n_h(h_values)

#get the trained n(h)
NN_n_h_values = NN_n_h(h_values, settings)

#@show size(h_values)
#@show h_values
#@show size(truth_n_h_values)
#@show truth_n_h_values
#@show size(NN_n_h_values)
#@show NN_n_h_values

#plot the results
plot_and_compare_n_h(h_values, truth_n_h_values, NN_n_h_values)

println("Done")
