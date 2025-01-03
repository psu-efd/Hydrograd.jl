
using Revise

ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0  # Disable auto precompilation

using Hydrograd

#change the working directory to the current directory
cd(dirname(@__FILE__))


function main()
    # Path to your SRH-2D files
    #srhhydro_file = "oneD_channel_with_bump_coarse.srhhydro"
    srhhydro_file = "simple.srhhydro"
    
    # Create SRH_2D_Data object
    srh_data = SRH_2D_Case.SRH_2D_Data(srhhydro_file)

    #dump the srh_data object
    dump(srh_data)
    
    # Get basic information about the case
    case_name = SRH_2D_Case.srh_2D_Data_get_case_name(srh_data)
    println("Case name: ", case_name)
    
    # Get mesh information
    num_elements, num_nodes = SRH_2D_Case.srh_2D_Data_get_mesh_size(srh_data)
    println("Mesh has $num_elements elements and $num_nodes nodes")
    
    # Get simulation time information
    start_time, end_time = SRH_2D_Case.srh_2D_Data_get_simulation_times(srh_data)
    dt = SRH_2D_Case.srh_2D_Data_get_simulation_time_step(srh_data)
    println("Simulation runs from $start_time to $end_time hours with dt = $dt seconds")

    # Get Manning's n values
    ManningN_cell = srh_data.ManningN_cell
    ManningN_node = srh_data.ManningN_node

    println("Manning's n values at cell centers: ", ManningN_cell)
    println("Manning's n values at nodes: ", ManningN_node)

    println("Done!")

end

# Run the example
main()