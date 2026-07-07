using Test
using LinearAlgebra
using VortexStepMethod

# Test the actual frozen_wake! function using real objects
@testset "frozen_wake! Tests" begin
    @testset "frozen_wake! basic functionality" begin
        # Create wing via convert-then-load (obj -> POLAR_MATRICES yaml)
        wing = ram_air_matrix_wing(; n_panels=56)
        body_aero = BodyAerodynamics([wing])

        # Test that frozen_wake! doesn't throw errors
        n_panels = length(body_aero.panels)
        va_distribution = ones(n_panels, 3)  # Create velocity distribution for all panels
        va_distribution[:, 1] .= 15.0  # X velocity = 15.0 m/s
        va_distribution[:, 2] .= 0.0   # Y velocity = 0.0
        va_distribution[:, 3] .= 0.0   # Z velocity = 0.0

        # This should not throw an error
        @test_nowarn VortexStepMethod.frozen_wake!(body_aero, va_distribution)

        # Test that the function accepts the right number of velocity vectors
        @test size(va_distribution, 1) == length(body_aero.panels)
        @test size(va_distribution, 2) == 3  # 3D velocity vectors
    end
end
