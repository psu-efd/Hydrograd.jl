#misc tools

function update_1d_array(arr, idx, values)
    #update a 1D array with a sparse update and avoid breaking the computational graph by creating a new array
    #arr is the array to update
    #idx is the array of indices to update
    #values is the array of values to update the array with    

    Zygote.ignore() do
        #println(" ")
        #println("update_1d_array")
        #println("Types:")
        #println("arr: ", typeof(arr))
        #println("idx: ", typeof(idx))
        #println("values: ", typeof(values))
    end
    
    # Create output array
    updated_arr = copy(arr)
    
    # Use broadcasting for the update
    updated_arr = [i ∈ idx ? values[findfirst(==(i), idx)] : arr[i] for i in 1:length(arr)]

    # Use non-mutating array comprehension
    #updated_arr = map(1:length(arr)) do i
    #    # Find if current index is in idx
    #    match_idx = findfirst(==(i), idx)
    #    # If found, use the new value; otherwise keep original
    #    isnothing(match_idx) ? arr[i] : values[match_idx]
    #end

    return updated_arr
end



