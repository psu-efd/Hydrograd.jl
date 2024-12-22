#misc tools

using Zygote

function update_1d_array(arr, idx, values)
    #update a 1D array with a sparse update and avoid breaking the computational graph by creating a new array
    #arr is the array to update
    #idx is the array of indices to update
    #values is the array of values to update the array with    
    
    # Create a new array with the updated value
    updated_arr = [i in idx ? values[findfirst(isequal(i), idx)] : arr[i] for i in 1:length(arr)]

    return updated_arr
end