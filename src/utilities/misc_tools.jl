#misc tools

function update_1d_array(arr, idx, values)
    #update a 1D array with a sparse update and avoid breaking the computational graph by creating a new array
    #arr is the array to update
    #idx is the array of indices to update
    #values is the array of values to update the array with    
    
    # Create a new array with the updated value
    #updated_arr = [i in idx ? values[findfirst(isequal(i), idx)] : arr[i] for i in 1:length(arr)]

    # Zygote.ignore() do
    #     println("arr = ", arr)
    #     println("idx = ", idx)
    #     println("values = ", values)
    # end

    # Create a dictionary for sparse updates
    #update_map = Dict(idx[i] => values[i] for i in 1:length(idx))

    # Generate a new array without modifying the original
    #updated_arr = [get(update_map, i, arr[i]) for i in 1:length(arr)]

    # Use non-mutating array comprehension
    updated_arr = map(1:length(arr)) do i
        # Find if current index is in idx
        match_idx = findfirst(==(i), idx)
        # If found, use the new value; otherwise keep original
        isnothing(match_idx) ? arr[i] : values[match_idx]
    end

    return updated_arr
end



