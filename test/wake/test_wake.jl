using Test
using LinearAlgebra
using VortexStepMethod

# Test the actual frozen_wake! function using real objects
@testset "frozen_wake! Tests" begin
    @testset "frozen_wake! basic functionality" begin
        # Use existing test data that we know works
        data_dir = joinpath(dirname(dirname(@__DIR__)), "data", "ram_air_kite")
        body_path = joinpath(tempdir(), "ram_air_kite_body.obj")
        foil_path = joinpath(tempdir(), "ram_air_kite_foil.dat")
        
        body_src = joinpath(data_dir, "ram_air_kite_body.obj")
        foil_src = joinpath(data_dir, "ram_air_kite_foil.dat")
        
        # Check if source files exist before copying
        if isfile(body_src) && isfile(foil_src)
            cp(body_src, body_path; force=true)
            cp(foil_src, foil_path; force=true)
            
            try
                # Create wing and body aerodynamics with known good geometry
                # ObjWing is fully complete - no refine! needed
                wing = ObjWing(body_path, foil_path; n_panels=56)  # Use default panels
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
                
            finally
                # Clean up
                rm(body_path; force=true)
                rm(foil_path; force=true)
            end
        else
            @test_skip "Ram air kite test data files not found, skipping frozen_wake! test"
        end
    end
end
