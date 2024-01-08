


# Assuming the SRHGEOM file format is known and consistent
struct RiverGeometry
    coordinates::Array{Float64,2}  # 2D array for coordinates (x, y)
    elevations::Vector{Float64}    # Array for elevations
    # Add other fields as necessary
end

"""
    readSRHGEOM(filename::String) -> RiverGeometry

Read a SRHGEOM file and return the geometric data.
"""
function readSRHGEOM(filename::String)
    # Check if file exists
    if !isfile(filename)
        error("File not found: $filename")
    end

    # Initialize data structures
    coords = []
    elevs = []

    # Open and process the file
    open(filename, "r") do file
        for line in eachline(file)
            # Assuming each line contains x, y coordinates and elevation, separated by spaces
            data = split(line)
            push!(coords, [parse(Float64, data[1]), parse(Float64, data[2])])
            push!(elevs, parse(Float64, data[3]))
        end
    end

    return RiverGeometry(coords, elevs)
end

