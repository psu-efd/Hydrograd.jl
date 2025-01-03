# Moduel to parse SRH-2D case

module SRH_2D_Case

using Revise

using DelimitedFiles

include("SRH_2D_SRHHydro.jl")
include("SRH_2D_SRHGeom.jl")
include("SRH_2D_SRHMat.jl")
include("SRH_2D_Data.jl")

# SRH_2D_SRHHydro
export SRH_2D_SRHHydro,
       srhhydro_get_case_name,
       srhhydro_get_hydroMat_fileName, 
       srhhydro_get_grid_file_name, 
       srhhydro_get_mat_file_name

# SRH_2D_SRHGeom
export SRH_2D_SRHGeom,
       srhgeom_compute_num_of_elements_nodes,
       srhgeom_get_num_of_elements_nodes,
       srhgeom_read_srhgeom_file,
       srhgeom_build_edges_and_boundary_edges,
       srhgeom_build_node_elements,
       srhgeom_get_vtk_cell_type_code

# SRH_2D_SRHMat
export SRH_2D_SRHMat,
       srhmat_find_cell_material_ID,
       srhmat_build_material_zones_data,
       srhmat_get_manningN_dict

# SRH_2D_Data
export SRH_2D_Data,
       srh_2D_Data_get_case_name,
       srh_2D_Data_get_mesh_size,
       srh_2D_Data_get_node_coordinates,
       srh_2D_Data_get_element_connectivity,
       srh_2D_Data_get_simulation_times,
       srh_2D_Data_get_simulation_time_step,
       srh_2D_Data_get_mesh_size,
       srh_2D_Data_get_node_coordinates,
       srh_2D_Data_get_element_connectivity


end # end of module SRH_2D_Case

