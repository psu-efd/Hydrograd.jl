using Zygote

# Test function that mimics your SWE structure
function test_reshape_vcat(parameters)

    x = [1.0, 2.0*parameters[1], 3.0, 4.0*parameters[2], 5.0]

    # Create a vector of vectors (like your updates)
    updates = [[x[i], x[i]^2, x[i]^3] for i in eachindex(x)]

    intermediate = reduce(vcat, updates)

    # Convert to matrix using reshape and vcat
#    result = reshape(intermediate, 5, :)

    result = [[intermediate[i], intermediate[i+1], intermediate[i+2]] for i in 1:3:length(intermediate)]


    Zygote.ignore() do
        @show typeof(updates)
        @show size(updates)
        @show updates
        @show typeof(intermediate)
        @show size(intermediate)
        @show intermediate
        @show typeof(result)
        @show size(result)
        @show result
    end
    
    # Return sum for testing gradient
    return sum(sum(x) for x in result)
end

# Alternative version using array comprehension for comparison
function test_comprehension(parameters)

    updates = rand(5,3).*parameters[1] .+ parameters[2]

    Zygote.ignore() do
        @show typeof(updates)
        @show size(updates)
        @show updates
    end
    

    result = [updates[i][j] for i in eachindex(updates), j in eachindex(updates[1])]

    Zygote.ignore() do
        @show typeof(result)
        @show size(result)
        @show result
    end
    
    return sum(sum(x) for x in result)
end

# Test data
parameters = [1.0, 2.0]

# Test forward pass
println("Testing reshape(reduce(vcat, ...))")
result1 = test_reshape_vcat(parameters)
println("Forward pass result: ", result1)
println("  ")

# Test gradient
grad1 = Zygote.gradient(test_reshape_vcat, parameters)[1]
println("Gradient 1: ", grad1)

# Compare with array comprehension version
println("\nTesting array comprehension")
result2 = test_comprehension(parameters)
println("Forward pass result: ", result2)
grad2 = Zygote.gradient(test_comprehension, parameters)[1]
println("Gradient 2: ", grad2)

# Verify gradients match
println("\nGradients match: ", all(isapprox.(grad1, grad2)))

# Verify with finite differences
function finite_diff(f, x, h=1e-7)
    grad = similar(x)
    for i in eachindex(x)
        x_plus = copy(x)
        x_plus[i] += h
        grad[i] = (f(x_plus) - f(x)) / h
    end
    return grad
end

fd_grad = finite_diff(test_reshape_vcat, parameters)
println("\nFinite difference gradient: ", fd_grad)
println("Gradient matches finite diff: ", all(isapprox.(grad1, fd_grad, rtol=1e-5)))
