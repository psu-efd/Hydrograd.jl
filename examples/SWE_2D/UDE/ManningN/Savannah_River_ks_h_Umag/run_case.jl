using Hydrograd

# Get the control file path
control_file = joinpath(@__DIR__, "run_control.json")

# Call the main solver function
elapsed_time = Hydrograd.solve_swe_2D(control_file)

println("Elapsed time: $elapsed_time seconds")

println("Done!")