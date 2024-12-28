using Zygote

#A notebook on arrays and AD

#1. 1D array and AD 

# Behavior of b = a for a 1D Array in Julia

#     If a is a 1D array, then b = a creates a reference to a. This means:
    
#         b points to the same memory as a.
#         Changes made to b will also affect a, and vice versa.

#         which is different from b = copy(a)

#         b = copy(a) creates a deep copy of a, so b and a are different objects in memory.
#         Changes made to b will not affect a, and vice versa.

println("1D array and AD")

# define a 1D array
a = [1.0, 2.0, 3.0]
@show a

# create a reference to a
b = a
@show b

# create a deep copy of a
c = copy(a)
@show c

# modify b
b[1] = 10.0
@show b
@show a

@show c

function f_1d_array_ref(arr)
    arr_ref = arr      # assign
    #arr_ref[1] = 10.0   #can't do array mutation; breaks AD
    return sum(arr_ref .^ 2) # Perform operations on the copy
end

function f_1d_array_copy(arr)
    arr_copy = copy(arr)      # assign
    arr_copy = arr_copy .^ 2
    return sum(arr_copy .^ 2) # Perform operations on the copy
end

grad_ref = Zygote.gradient(f_1d_array_ref, [1.0, 2.0, 3.0])  # Output: [2.0, 4.0, 6.0] (works fine)
grad_copy = Zygote.gradient(f_1d_array_copy, [1.0, 2.0, 3.0])  # Output: [2.0, 4.0, 6.0] (works fine)

println("grad_ref: $grad_ref")
println("grad_copy: $grad_copy")

#2. 2D array and AD
println("2D array and AD")

# define a 2D array
a = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
@show a

# create a reference to a
b = a
@show b

# create a deep copy of a
c = copy(a)
@show c 

# create a view of a
d = view(a, :, 1)
@show d

# modify d
d[1] = 10.0
@show d
@show a

# create a slice of 2D array creates a new array
println("create a slice of a")
a = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
@show a
e = a[:, 1]
@show e
e[1] = 10.0 # this modifies e, but not a
@show e
@show a  # a is not modified

function f_2d_array_ref(arr)
    arr_ref = arr      # assign
    #arr_ref[1, 1] = 10.0   #can't do array mutation; breaks AD
    return sum(arr_ref .^ 2) # Perform operations on the copy
end

function f_2d_array_copy(arr)   
    arr_copy = copy(arr)      # assign
    arr_copy = arr_copy .^ 2
    return sum(arr_copy .^ 2) # Perform operations on the copy
end

function f_2d_array_view(arr)
    arr_view = view(arr, :, 1)
    #arr_view[1] = 10.0
    return sum(arr_view .^ 2) # Perform operations on the copy
end

function f_2d_array_slice(arr)
    arr_slice = arr[:, 1]
    arr_slice = arr_slice .* 2.0
    return sum(arr_slice .^ 2) # Perform operations on the copy
end

param = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
@show param

grad_ref_2d = Zygote.gradient(f_2d_array_ref, param)  
grad_copy_2d = Zygote.gradient(f_2d_array_copy, param) 
grad_view_2d = Zygote.gradient(f_2d_array_view, param) 
grad_slice_2d = Zygote.gradient(f_2d_array_slice, param) 

println("grad_ref_2d: $grad_ref_2d")
println("grad_copy_2d: $grad_copy_2d")
println("grad_view_2d: $grad_view_2d")
println("grad_slice_2d: $grad_slice_2d")

#test whether I can modify an argument in a function to be differentiated. 

misc_arg = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
arr = [1.0, 2.0, 3.0]

function f_modify_arg(arr, misc_arg)
    arr_ref = arr      # assign
    misc_arg = misc_arg .* 2.0
    return sum(arr_ref .^ 2) # Perform operations on the copy
end

grad_modify_arg = Zygote.gradient((x) -> f_modify_arg(x, misc_arg), arr)

println("grad_modify_arg: $grad_modify_arg")
println("misc_arg: $misc_arg")










