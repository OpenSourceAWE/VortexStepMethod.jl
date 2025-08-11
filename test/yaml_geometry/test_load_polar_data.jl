using VortexStepMethod
using VortexStepMethod: load_polar_data
using LinearAlgebra
using Test
using YAML
using Logging

@testset "load_polar_data Function Tests" begin
    # Setup temporary files for testing
    test_csv_path = joinpath(tempdir(), "test_polar.csv")
    
    # Clean up function
    function cleanup_test_files()
        isfile(test_csv_path) && rm(test_csv_path; force=true)
    end
    
    @testset "Valid CSV File" begin
        # Create a valid CSV file with polar data
        csv_content = """alpha,cl,cd,cm
              -10.0,0.1,0.02,-0.01
              -5.0,0.5,0.015,-0.005
              0.0,1.0,0.01,0.0
              5.0,1.5,0.015,0.005
              10.0,1.8,0.025,0.01
              15.0,1.6,0.05,0.015"""
        write(test_csv_path, csv_content)
        
        aero_data, aero_model = load_polar_data(test_csv_path)
        
        @test aero_model == POLAR_VECTORS
        @test aero_data isa Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
        @test length(aero_data) == 4  # alpha, cl, cd, cm
        @test length(aero_data[1]) == 6  # 6 data points
        
        # Test that alpha is converted to radians
        @test aero_data[1] ≈ deg2rad.([-10.0, -5.0, 0.0, 5.0, 10.0, 15.0])
        
        # Test coefficient values
        @test aero_data[2] ≈ [0.1, 0.5, 1.0, 1.5, 1.8, 1.6]  # cl
        @test aero_data[3] ≈ [0.02, 0.015, 0.01, 0.015, 0.025, 0.05]  # cd
        @test aero_data[4] ≈ [-0.01, -0.005, 0.0, 0.005, 0.01, 0.015]  # cm
    end
    
    @testset "Different Column Order" begin
        # Test with different column ordering
        csv_content = """cm,alpha,cd,cl
                -0.01,-10.0,0.02,0.1
                -0.005,-5.0,0.015,0.5
                0.0,0.0,0.01,1.0"""
        write(test_csv_path, csv_content)
        
        aero_data, aero_model = load_polar_data(test_csv_path)
        
        @test aero_model == POLAR_VECTORS
        @test aero_data[1] ≈ deg2rad.([-10.0, -5.0, 0.0])  # alpha
        @test aero_data[2] ≈ [0.1, 0.5, 1.0]  # cl
        @test aero_data[3] ≈ [0.02, 0.015, 0.01]  # cd
        @test aero_data[4] ≈ [-0.01, -0.005, 0.0]  # cm
    end
    
    @testset "Case Insensitive Headers" begin
        # Test with mixed case headers
        csv_content = """ALPHA,CL,CD,CM
          0.0,1.0,0.01,0.0
          5.0,1.5,0.015,0.005"""
        write(test_csv_path, csv_content)
        
        aero_data, aero_model = load_polar_data(test_csv_path)
        
        @test aero_model == POLAR_VECTORS
        @test length(aero_data[1]) == 2
    end
    
    @testset "Missing File" begin
        nonexistent_path = joinpath(tempdir(), "nonexistent.csv")
        aero_data, aero_model = suppress_warnings(() -> load_polar_data(nonexistent_path))
        
        @test aero_model == INVISCID
        @test aero_data === nothing
    end
    
    @testset "Empty File Path" begin
        aero_data, aero_model = suppress_warnings(() -> load_polar_data(""))
        
        @test aero_model == INVISCID
        @test aero_data === nothing
    end
    
    @testset "Invalid CSV Format" begin
        # Missing required columns
        csv_content = """alpha,cl
        0.0,1.0
        5.0,1.5"""
        write(test_csv_path, csv_content)
        
        aero_data, aero_model = suppress_warnings(() -> load_polar_data(test_csv_path))
        
        @test aero_model == INVISCID
        @test aero_data === nothing
    end
    
    @testset "Empty CSV File" begin
        write(test_csv_path, "")
        
        aero_data, aero_model = suppress_warnings(() -> load_polar_data(test_csv_path))
        
        @test aero_model == INVISCID
        @test aero_data === nothing
    end
    
    @testset "Malformed Data" begin
        # Invalid numeric data
        csv_content = """alpha,cl,cd,cm
        invalid,1.0,0.01,0.0
        5.0,1.5,0.015,0.005"""
        write(test_csv_path, csv_content)
        
        aero_data, aero_model = suppress_warnings(() -> load_polar_data(test_csv_path))
        
        @test aero_model == INVISCID
        @test aero_data === nothing
    end
    
    # Cleanup after all tests
    cleanup_test_files()
end
