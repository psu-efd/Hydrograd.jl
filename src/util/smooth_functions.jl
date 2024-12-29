#smooth functions: To make the numerical derivatives of these functions well defined and/or more stale

#smooth the absolute value function
function smooth_abs(x)
    return sqrt(x^2 + eps(typeof(x)))
end

#smooth the sign function
function smooth_sign(x)
    return x / smooth_abs(x)
end

#smooth the max function
function smooth_max(x, y)
    return (x + y + smooth_abs(x - y)) / 2.0
end

#smooth the min function
function smooth_min(x, y)
    return (x + y - smooth_abs(x - y)) / 2.0
end

#smooth the sqrt function
function smooth_sqrt(x)
    return sqrt(x + eps(typeof(x)))
end

#smooth the pow function
function smooth_pow(x, y)
    return (x + eps(typeof(x)))^y  
end


