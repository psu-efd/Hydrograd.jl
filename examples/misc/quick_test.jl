using Zygote

# Your case:
#my_function(u, p)  # u is q_x (12), p is ManningN (12)
# Output: 12×3 matrix
# Gradient wrt p: vector of 12 (matches p's dimension)

# Simple example to illustrate:
function example(p)
    # p is a vector of length 3
    # Returns a 3×2 matrix
    return [p[i] * j for i in 1:3, j in 1:2]
end

x = [1.0, 2.0, 3.0]
_dy, back = Zygote.pullback(example, x)
@show size(_dy)    # (3, 2) - function output
λ = ones(3, 2)     # matches function output
grad = back(λ)[1]
@show size(grad)   # (3,) - matches input parameter size