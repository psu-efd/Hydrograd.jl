using Revise

using AdHydraulics

# Example usage
# filename = joinpath(@__DIR__, "minimal_for_debug.msh")
# nodes, elements = read_gmsh_v2(filename)
# element_areas, element_centroids, element_normals = compute_mesh_properties(nodes, elements)

# # Output results
# println("Nodes:")
# println(nodes)
# println("Elements:")
# println(elements)
# println("Element Areas:")
# println(element_areas)
# println("Element Centroids:")
# println(element_centroids)
# println("Element Normals:")
# println(element_normals)

# Example usage
# file_path = joinpath(@__DIR__, "backwater.srhgeom" ) 
# name, grid_unit, nodes, elements, boundaries = read_srhgeom(file_path)

# # Compute geometric and topological properties
# element_areas, element_centroids, element_normals = compute_mesh_properties_srhgeom(nodes, elements)

# # Output results
# println("Mesh Name: $name")
# println("Grid Unit: $grid_unit")
# println("Number of Nodes: $(length(nodes))")
# println("Number of Elements: $(length(elements))")
# println("Number of Boundaries: $(length(boundaries))")
# println("boundaries: ", boundaries)
# println("First 5 Elements: ", Dict(first(elements, 5)))
# println("First 5 Node Coordinates: ", Dict(first(nodes, 5)))
# println("First 5 Element Areas: ", Dict(first(element_areas, 5)))
# println("First 5 Element Centroids: ", Dict(first(element_centroids, 5)))


# Example usage
# file_path = joinpath(@__DIR__, "backwater.srhmat" ) 
# num_materials, material_names, material_elements = read_srhmat(file_path)

# # Output parsed information
# println("Number of Materials: $num_materials")
# println("Material Names: ", material_names)
# println("Material Elements (First Material): ", first(material_elements, 1))


# Path to the srhhydro file
file_path = joinpath(@__DIR__, "backwater.srhhydro" ) 

# # Call the function and display the parsed data
# parsed_data = read_srhhydro(file_path)
# for (key, value) in parsed_data
#     println("$key => $value")
# end

#using .SRH2D  # Load the module (ensure it is in the current environment or included)

# Initialize the SRH_2D_Data object
srh_data = SRH2D.SRH_2D_Data(file_path)

# Access data from the srhhydro file
println("Case Name: ", srh_data.hydro.content["Case"])
println("Grid File: ", srh_data.hydro.content["Grid"])
println("HydroMat File: ", srh_data.hydro.content["HydroMat"])
println("Manning's N Values: ", srh_data.hydro.content["ManningsN"])

# Access the geom file information
println("Number of Nodes: ", srh_data.geom.num_nodes)
println("Number of Elements: ", srh_data.geom.num_elements)
println("Node Coordinates: ", srh_data.geom.nodes)
println("Element Connectivity: ", srh_data.geom.elements)

# Access material data
println("Material Zones: ", srh_data.mat.materials)

# Modify data in the srhhydro file
srh_data.hydro.content["Case"] = "New Case Name"
println("Updated Case Name: ", srh_data.hydro.content["Case"])

# Save modifications to a new srhhydro file
new_hydro_file = "path/to/updated_backwater.srhhydro"
open(new_hydro_file, "w") do io
    for (key, value) in srh_data.hydro.content
        if typeof(value) <: Dict
            for (subkey, subvalue) in value
                println(io, "$key $subkey $subvalue")
            end
        else
            println(io, "$key $value")
        end
    end
end

# Example query
println("Boundary Conditions: ", srh_data.hydro.content["BC"])
println("IQ Parameters: ", srh_data.hydro.content["IQParams"])
