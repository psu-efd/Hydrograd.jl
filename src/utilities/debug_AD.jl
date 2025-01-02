#code piece for debugging AD correctness (by directly inclue this file in the code)

function debug_AD(ode_f, Q0, swe_2D_constants, params_vector, swe_extra_params)

    #debug start
    @show Q0
    @show swe_2D_constants.tspan
    @show params_vector
    @assert !ismissing(params_vector) "params_vector contains missing values!"
    @assert !isnothing(params_vector) "params_vector contains `nothing` values!"

    prob = ODEProblem(ode_f, Q0, swe_2D_constants.tspan, params_vector)

    #use Enzyme to test the gradient of the ODE and identify the source of the error
    #See https://docs.sciml.ai/SciMLSensitivity/dev/faq/
    SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true
    p = prob.p
    y = prob.u0
    f = prob.f
    t = swe_2D_constants.tspan[1]  # Add this line to define t
    @show typeof(p)
    @show typeof(y)
    @show typeof(f)
    @show typeof(t)
    @show p
    @show y
    @show f
    @show t

    # Test forward pass first
    try
        #test_forward = f(y, p, t, swe_extra_params)
        
        #test_forward = swe_2d_rhs(y, p, t, swe_extra_params)

        #@show typeof(test_forward)
        #@show size(test_forward)
        #@show test_forward

        println("\nForward pass successful\n")
        
    catch e
        println("\nForward pass failed\n")
        @show e

        #stop here
        @assert false "stop after forward pass debug"
    end

    # Now test the pullback with more detailed error catching
    try
        λ = ones(size(prob.u0)) #zero(prob.u0)

        #Zygote.pullback takes two arguments:
        #  First argument: a function that we want to differentiate
        #  Remaining arguments: the values at which to evaluate the function (y and p in this case)
        #_dy is the result of the forward pass; back is the gradient function
        #_dy, back = Zygote.pullback((u, p) -> f(u, p, t), y, p)
        #_dy, back = Zygote.pullback((u, p) -> Array(swe_2d_rhs(u, p, t, swe_extra_params)), y, p)
        _dy, back = Zygote.pullback((u, p) -> swe_2d_rhs(u, p, t, swe_extra_params), y, p)

        #_dy, back = Zygote.pullback(y, p) do u, p  
        #    #vec(f(u, p, t))
        #    f(u, p, t)
        #end
        println("\nPullback creation successful\n")
        @show typeof(_dy)
        @show size(_dy)
        @show _dy

        try
             # Convert λ to match _dy type
            λ = convert(typeof(_dy), λ)

            tmp1, tmp2 = back(λ)                  #tmp1 is the gradient of the state variables; tmp2 is the gradient of the parameters
            println("\nBackward pass is successful\n")
            @show typeof(tmp1)
            @show size(tmp1)
            @show typeof(tmp2)
            @show size(tmp2)
            @show tmp1
            @show tmp2
        catch e
            println("\nBackward pass failed\n")
            @show e
            @show typeof(λ)
            @show size(λ)
            @show λ

            #stop here
            #return
            #@assert false "stop after backward pass debug"
        end
    catch e
        println("\nPullback creation failed\n")
        @show e

        #stop here
        @assert false "stop after pullback creation debug"
    end

    #stop here
    @assert false "stop after AD debug"

end



