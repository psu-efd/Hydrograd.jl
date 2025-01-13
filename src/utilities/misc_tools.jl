#misc tools

# This version make a copy of values, which is the full lenght of the array; effectively to reorder the values based on the indices
function update_1d_array(idx::AbstractVector{Int}, values::AbstractVector{T1}) where {T1}
    #update a 1D array with a sparse update and avoid breaking the computational graph by creating a new array
    #arr is the array to update
    #idx is the array of indices to update
    #values is the array of values to update the array with    
   
    # Create output array
    #updated_arr = copy(values)
    
    # Use broadcasting for the update
    #updated_arr = [values[idx[i]] for i in eachindex(idx)]
    updated_arr = [values[idx[i]] for i in eachindex(idx)]

    return updated_arr
end

# function update_1d_array(arr::AbstractVector{T1}, idx::AbstractVector{Int}, values::AbstractVector{T2}) where {T1, T2}
#     #update a 1D array with a sparse update and avoid breaking the computational graph by creating a new array
#     #arr is the array to update
#     #idx is the array of indices to update
#     #values is the array of values to update the array with    
   
#     # Create output array
#     updated_arr = copy(arr)
    
#     # Use broadcasting for the update
#     updated_arr = [i ∈ idx ? values[findfirst(==(i), idx)] : arr[i] for i in eachindex(arr)]

#     return updated_arr
# end

# When saving the ODE solution, only save the solution data
function save_ode_solution(sol, filename)
    # Extract just the solution data
    sol_data = Dict(
        "t" => sol.t,
        "u" => sol.u,
        "retcode" => sol.retcode,
        "stats" => Dict(
            "nf" => sol.stats.nf,
            "naccept" => sol.stats.naccept,
            "nreject" => sol.stats.nreject
        )
    )
    
    # Save only the solution data
    jldsave(filename; sol_data=sol_data)
end





