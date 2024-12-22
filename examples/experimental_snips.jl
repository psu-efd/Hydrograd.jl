using Zygote

function update_array(arr, idx, value)
    # Create a new array with the updated value
    updated_arr = arr .* (1 .- idx) .+ value .* idx
    return updated_arr
end

# Example: Update the second element of the array
arr = [1.0, 2.0, 3.0]
idx = [0, 1, 0]  # Indicator for which element to update
value = 10.0

updated_arr = update_array(arr, idx, value)
println(updated_arr)  # Output: [1.0, 10.0, 3.0]

# Works with Zygote.gradient
grad = Zygote.gradient(x -> sum(update_array(x, idx, value)), arr)
println(grad)
