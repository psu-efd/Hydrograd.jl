using SciMLBase
using OrdinaryDiffEq
using Zygote, ForwardDiff

# Test different solvers
println("Tsit5: ", SciMLBase.isautodifferentiable(Tsit5()))
println("Euler: ", SciMLBase.isautodifferentiable(Euler()))

# Test different sensitivity algorithms
#println("ZygoteVJP: ", SciMLBase.isautodifferentiable(ZygoteVJP()))
#println("ForwardDiffSensitivity: ", SciMLBase.isautodifferentiable(ForwardDiffSensitivity()))
#println("InterpolatingAdjoint: ", SciMLBase.isautodifferentiable(InterpolatingAdjoint()))