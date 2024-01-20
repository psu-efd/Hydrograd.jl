using AdHydraulics

using DifferentialEquations, Optimization, OptimizationPolyalgorithms,
    Zygote, SciMLSensitivity

using ForwardDiff


"""
function loss(θ)
    return sum(abs.(θ))
end 

ps = [1.2, -8.4, 2.5, 0.0]

grad = ForwardDiff.gradient(loss, ps)
"""

L = 1.8
nCells = 5
mesh = initialize_mesh_1D(nCells, L)

println(mesh)

swe_1d_para = swe_1D_parameters(
    g=9.81,
    ManningN=0.03,
    t=0.0,
    dt_min=0.005,
    dt=0.005,
    CFL=0.001,
    tEnd=0.0051,
    h_small=1.0e-3,
    RiemannSolver="HLL"
)

#  create the left and right boundaries
leftBoundary = Boundary_1D(bcType=inletQ::Boundary_Type_Name, bcValue=0.5)
rigthBoundary = Boundary_1D(bcType=exitH::Boundary_Type_Name, bcValue=1.0)

#leftBoundary = Boundary_1D(bcType=wall::Boundary_Type_Name, bcValue=0.0)
#rigthBoundary = Boundary_1D(bcType=wall::Boundary_Type_Name, bcValue=0.0)


fields = initialize_swe_1D_fields(mesh, swe_1d_para, leftBoundary, rigthBoundary)

println(fields)