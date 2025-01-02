using Infiltrator

function compute_area(radius)
    area = π * radius^2
    @infiltrate # Pause execution here
    return area
end

result = compute_area(3)
println("Computed area: $result")
