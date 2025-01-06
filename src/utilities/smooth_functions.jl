#smooth functions: To make the numerical derivatives of these functions well defined and/or more stale

#smooth the absolute value function
function smooth_abs(x::T) where T <: Real
    # Check if x is a real number
    #@assert eltype(x) <: Real "x must be a real number"

    # Ensure x is not an array
    #@assert !isa(x, AbstractArray) "x must not be an array"

    return sqrt(x^2 + eps(T))
end

#smooth the sign function
function smooth_sign(x::T) where T <: Real
    # Ensure x is not an array
    #@assert !isa(x, AbstractArray) "x must not be an array"

    return x / smooth_abs(x)
end

#smooth the max function
function smooth_max(x::T, y::T) where T <: Real
    # Ensure x is not an array
    #@assert !isa(x, AbstractArray) "x must not be an array"

    return (x + y + smooth_abs(x - y)) / 2.0
end

#smooth the min function
function smooth_min(x::T, y::T) where T <: Real
    # Ensure x is not an array
    #@assert !isa(x, AbstractArray) "x must not be an array"

    return (x + y - smooth_abs(x - y)) / 2.0
end

#smooth the sqrt function
function smooth_sqrt(x::T) where T <: Real
    # Ensure x is not an array
    #@assert !isa(x, AbstractArray) "x must not be an array"

    return sqrt(x + eps(T))
end

#smooth the pow function
function smooth_pow(x::T, y::Int) where T <: Real
    # Ensure x is not an array
    #@assert !isa(x, AbstractArray) "x must not be an array"
    
    return (x + eps(T))^y  
end


