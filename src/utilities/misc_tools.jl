#misc tools

function update_1d_array(arr::AbstractVector{T}, idx::AbstractVector{Int}, values::AbstractVector{T}) where T<:Real
    #update a 1D array with a sparse update and avoid breaking the computational graph by creating a new array
    #arr is the array to update
    #idx is the array of indices to update
    #values is the array of values to update the array with    
   
    # Create output array
    updated_arr = copy(arr)
    
    # Use broadcasting for the update
    updated_arr = [i ∈ idx ? values[findfirst(==(i), idx)] : arr[i] for i in eachindex(arr)]

    return updated_arr
end

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





