using JSON3

# Define the file paths
save_path = dirname(@__FILE__)
json_file = joinpath(save_path, "inversion_config.json")
output_file = joinpath(save_path, "modified_config.json")


# Read the JSON file
config = JSON3.read(open(json_file), Dict)

# Print the original configuration
println("Original Configuration:")
println(config)

# Modify the configuration
# Example modifications:
# 1. Change the learning rate in the "model" section
config["model"]["learning_rate"] = 0.01

# 2. Change the batch size in the "trainer" section
config["trainer"]["batch_size"] = 20

# 3. Add a new field to "inverter"
config["inverter"]["new_parameter"] = "new_value"

# Print the modified configuration
println("Modified Configuration:")
println(config)

# Save the modified JSON back to a file
open(output_file, "w") do io
    JSON3.pretty(io, config)
    println(io)
    #JSON3.write(io, config; indent=4, pretty=true)
end

println("Modified configuration saved to $output_file.")
