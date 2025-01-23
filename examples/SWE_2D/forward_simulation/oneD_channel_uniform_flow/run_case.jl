using Hydrograd

#change the current working directory to the directory of the file
current_dir = dirname(@__FILE__)
cd(current_dir)

# Get the control file path
control_file = joinpath(@__DIR__, "run_control.json")

# Call the main solver function
elapsed_time = Hydrograd.solve_swe_2D(control_file)

println("Elapsed time: $elapsed_time seconds")

#change the current working directory back to the original directory
cd(current_dir)

println("Done!")