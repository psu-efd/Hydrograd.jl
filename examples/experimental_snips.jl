using Zygote

list_of_vectors = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]  # List of vectors

# Create a matrix from the list of vectors
matrix = [list_of_vectors[i][j] for i in 1:length(list_of_vectors), j in 1:length(list_of_vectors[1])]

display(matrix)

# Test with Zygote gradient
grad = Zygote.gradient(x -> sum([x[i][j] for i in 1:length(x), j in 1:length(x[1])]), list_of_vectors)

println(grad)


