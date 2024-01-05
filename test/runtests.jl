using AdHydraulics
using Test

@testset "AdHydraulics.jl" begin
    # Write your tests here.
    @test AdHydraulics.greet_your_package_name() == "Hello from AdHydraulics!"
    @test AdHydraulics.greet_your_package_name() != "Hello world!"
end
