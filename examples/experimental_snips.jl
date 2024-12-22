list_of_vectors = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]  # List of vectors

array = vcat(list_of_vectors'...)  # Use transpose (') for proper alignment
println(array)
display(array)
println(typeof(array))
