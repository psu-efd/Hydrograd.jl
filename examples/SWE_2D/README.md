
# Running Instructions

1. Open the Julia REPL, activate the Hydrograd environment, and navigate to the case directory, e.g.

   ```bash
   julia> cd("examples/SWE_2D/forward_simulation/oneD_channel_with_bump")
   ```

2. Run the script by typing
   ```bash
   julia> include("run_case.jl")
   ```
   The script `run_case.jl` is located in each of the case directories. This script is pretty much uniform for all cases and all applications. It reades the `run_control.json` file within the case directory to get the parameters for the case run, and then calls the `solve_swe_2D.jl` script to solve the shallow water equations and perform other tasks.
3. The script will execute and generate the output files in the current case directory.
4. Post-processing can be done using the Python scripts in the case directory. Most of the cases also output the results in VTK format, which can be visualized using [ParaView](https://www.paraview.org/).

