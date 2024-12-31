# Create a RefValue
x = Base.RefValue{Any}(10)  # or Base.RefValue{Any}(10)

# Access the value using []
println(x[])  # prints: 10

# Modify the value
x[] = "hello"
println(x[])  # prints: "hello"