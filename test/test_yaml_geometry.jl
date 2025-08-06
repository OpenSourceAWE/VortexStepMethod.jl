using VortexStepMethod
using VortexStepMethod: load_polar_data, YamlWing
using LinearAlgebra
using Test
using YAML

@testset "YAML Geometry Tests" begin
    # Setup temporary files for testing
    test_csv_path = joinpath(tempdir(), "test_polar.csv")
    test_yaml_path = joinpath(tempdir(), "test_wing.yaml")
    test_polar_dir = joinpath(tempdir(), "polars")
    
    # Clean up function
    function cleanup_test_files()
        for file in [test_csv_path, test_yaml_path]
            isfile(file) && rm(file)
        end
        isdir(test_polar_dir) && rm(test_polar_dir; recursive=true)
    end
    
    @testset "load_polar_data Function Tests" begin
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
            aero_data, aero_model = load_polar_data(nonexistent_path)
            
            @test aero_model == INVISCID
            @test aero_data === nothing
        end
        
        @testset "Empty File Path" begin
            aero_data, aero_model = load_polar_data("")
            
            @test aero_model == INVISCID
            @test aero_data === nothing
        end
        
        @testset "Invalid CSV Format" begin
            # Missing required columns
            csv_content = """alpha,cl
0.0,1.0
5.0,1.5"""
            write(test_csv_path, csv_content)
            
            aero_data, aero_model = load_polar_data(test_csv_path)
            
            @test aero_model == INVISCID
            @test aero_data === nothing
        end
        
        @testset "Empty CSV File" begin
            write(test_csv_path, "")
            
            aero_data, aero_model = load_polar_data(test_csv_path)
            
            @test aero_model == INVISCID
            @test aero_data === nothing
        end
        
        @testset "Malformed Data" begin
            # Invalid numeric data
            csv_content = """alpha,cl,cd,cm
invalid,1.0,0.01,0.0
5.0,1.5,0.015,0.005"""
            write(test_csv_path, csv_content)
            
            aero_data, aero_model = load_polar_data(test_csv_path)
            
            @test aero_model == INVISCID
            @test aero_data === nothing
        end
    end
    
    @testset "YamlWing Constructor Tests" begin
        # Create polar data directory and files
        mkpath(test_polar_dir)
        
        # Create polar files
        polar1_path = joinpath(test_polar_dir, "1.csv")
        polar2_path = joinpath(test_polar_dir, "2.csv")
        
        polar1_content = """alpha,cl,cd,cm
-10.0,0.1,0.02,-0.01
0.0,1.0,0.01,0.0
10.0,1.8,0.025,0.01"""
        
        polar2_content = """alpha,cl,cd,cm
-10.0,0.15,0.025,-0.015
0.0,1.1,0.012,0.001
10.0,1.9,0.03,0.012"""
        
        write(polar1_path, polar1_content)
        write(polar2_path, polar2_content)
        
        @testset "Valid YAML Wing Construction" begin
            # Create a valid YAML wing configuration
            yaml_content = """
wing_sections:
  headers: [airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]
  data:
    - [1, 0.0, 4.08, 8.05, 2.0, 4.08, 8.05]
    - [2, 0.0, -4.08, 8.05, 2.0, -4.08, 8.05]
    - [1, 0.0, 0.0, 4.9, 2.0, 0.0, 4.9]

wing_airfoils:
  alpha_range: [-10, 10, 1]
  reynolds: 1000000
  headers: [airfoil_id, type, info_dict]
  data:
    - [1, polars, {csv_file_path: "polars/1.csv"}]
    - [2, polars, {csv_file_path: "polars/2.csv"}]
"""
            write(test_yaml_path, yaml_content)
            
            wing = YamlWing(test_yaml_path; n_panels=4, n_groups=2)
            
            @test wing isa Wing
            @test wing.n_panels == 4
            @test wing.n_groups == 2
            @test wing.spanwise_distribution == LINEAR
            @test wing.spanwise_direction ≈ [0.0, 1.0, 0.0]
            @test length(wing.sections) == 3
            
            # Test section coordinates
            @test wing.sections[1].LE_point ≈ [0.0, 4.08, 8.05]
            @test wing.sections[1].TE_point ≈ [2.0, 4.08, 8.05]
            @test wing.sections[2].LE_point ≈ [0.0, -4.08, 8.05]
            @test wing.sections[2].TE_point ≈ [2.0, -4.08, 8.05]
            @test wing.sections[3].LE_point ≈ [0.0, 0.0, 4.9]
            @test wing.sections[3].TE_point ≈ [2.0, 0.0, 4.9]
            
            # Test that sections have polar data
            @test wing.sections[1].aero_model == POLAR_VECTORS
            @test wing.sections[2].aero_model == POLAR_VECTORS
            @test wing.sections[3].aero_model == POLAR_VECTORS
            
            # Test polar data is loaded
            @test wing.sections[1].aero_data isa Tuple
            @test length(wing.sections[1].aero_data) == 4
            @test length(wing.sections[1].aero_data[1]) == 3  # 3 alpha points
        end
        
        @testset "Wing Constructor Parameters" begin
            yaml_content = """
wing_sections:
  headers: [airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]
  data:
    - [1, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0]
    - [1, 0.0, -1.0, 0.0, 1.0, -1.0, 0.0]

wing_airfoils:
  alpha_range: [-10, 10, 1]
  reynolds: 1000000
  headers: [airfoil_id, type, info_dict]
  data:
    - [1, polars, {csv_file_path: "polars/1.csv"}]
"""
            write(test_yaml_path, yaml_content)
            
            # Test custom parameters
            wing = YamlWing(
                test_yaml_path;
                n_panels=8,
                n_groups=4,
                spanwise_distribution=COSINE,
                remove_nan=false
            )
            
            @test wing.n_panels == 8
            @test wing.n_groups == 4
            @test wing.spanwise_distribution == COSINE
            @test !wing.remove_nan
        end
        
        @testset "Missing Polar Files" begin
            # Create YAML with non-existent polar files
            yaml_content = """
wing_sections:
  headers: [airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]
  data:
    - [1, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0]

wing_airfoils:
  alpha_range: [-10, 10, 1]
  reynolds: 1000000
  headers: [airfoil_id, type, info_dict]
  data:
    - [1, polars, {csv_file_path: "nonexistent.csv"}]
"""
            write(test_yaml_path, yaml_content)
            
            wing = YamlWing(test_yaml_path; n_panels=2)
            
            # Should fall back to INVISCID model
            @test wing.sections[1].aero_model == INVISCID
            @test wing.sections[1].aero_data === nothing
        end
        
        @testset "Sections Without Polar Files" begin
            # Create YAML with empty polar file paths
            yaml_content = """
wing_sections:
  headers: [airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]
  data:
    - [1, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0]

wing_airfoils:
  alpha_range: [-10, 10, 1]
  reynolds: 1000000
  headers: [airfoil_id, type, info_dict]
  data:
    - [1, polars, {csv_file_path: ""}]
"""
            write(test_yaml_path, yaml_content)
            
            wing = YamlWing(test_yaml_path; n_panels=2)
            
            # Should fall back to INVISCID model
            @test wing.sections[1].aero_model == INVISCID
        end
        
        @testset "Invalid Parameters" begin
            yaml_content = """
wing_sections:
  headers: [airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]
  data:
    - [1, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0]

wing_airfoils:
  alpha_range: [-10, 10, 1]
  reynolds: 1000000
  headers: [airfoil_id, type, info_dict]
  data:
    - [1, polars, {csv_file_path: "polars/1.csv"}]
"""
            write(test_yaml_path, yaml_content)
            
            # Test invalid n_panels/n_groups combination
            @test_throws ArgumentError YamlWing(test_yaml_path; n_panels=5, n_groups=2)
            
            # Test invalid spanwise direction
            @test_throws ArgumentError YamlWing(test_yaml_path; spanwise_direction=[1.0, 0.0, 0.0])
        end
        
        @testset "Relative Path Resolution" begin
            # Test that relative paths in CSV files are resolved relative to YAML file
            subdir = joinpath(tempdir(), "subtest")
            mkpath(subdir)
            subdir_polar_dir = joinpath(subdir, "polars")
            mkpath(subdir_polar_dir)
            
            # Create polar file in subdirectory
            subdir_polar_path = joinpath(subdir_polar_dir, "airfoil.csv")
            write(subdir_polar_path, polar1_content)
            
            # Create YAML in subdirectory with relative path
            subdir_yaml_path = joinpath(subdir, "wing.yaml")
            yaml_content = """
wing_sections:
  headers: [airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]
  data:
    - [1, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0]

wing_airfoils:
  alpha_range: [-10, 10, 1]
  reynolds: 1000000
  headers: [airfoil_id, type, info_dict]
  data:
    - [1, polars, {csv_file_path: "polars/airfoil.csv"}]
"""
            write(subdir_yaml_path, yaml_content)
            
            wing = YamlWing(subdir_yaml_path; n_panels=2)
            
            # Should successfully load polar data with relative path
            @test wing.sections[1].aero_model == POLAR_VECTORS
            @test wing.sections[1].aero_data isa Tuple
            
            # Cleanup
            rm(subdir; recursive=true)
        end
        
        @testset "Complex Wing Geometry" begin
            # Test with more complex wing geometry
            yaml_content = """
wing_sections:
  headers: [airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]
  data:
    - [1, 0.0, 5.0, 0.0, 2.0, 5.0, 0.5]
    - [2, 0.2, 3.0, 0.2, 1.8, 3.0, 0.3]
    - [2, 0.4, 1.0, 0.4, 1.6, 1.0, 0.1]
    - [1, 0.5, 0.0, 0.5, 1.5, 0.0, 0.0]
    - [2, 0.4, -1.0, 0.4, 1.6, -1.0, 0.1]
    - [2, 0.2, -3.0, 0.2, 1.8, -3.0, 0.3]
    - [1, 0.0, -5.0, 0.0, 2.0, -5.0, 0.5]

wing_airfoils:
  alpha_range: [-10, 10, 1]
  reynolds: 1000000
  headers: [airfoil_id, type, info_dict]
  data:
    - [1, polars, {csv_file_path: "polars/1.csv"}]
    - [2, polars, {csv_file_path: "polars/2.csv"}]
"""
            write(test_yaml_path, yaml_content)
            
            wing = YamlWing(test_yaml_path; n_panels=12, n_groups=3)
            
            @test wing.n_panels == 12
            @test wing.n_groups == 3
            @test length(wing.sections) == 7
            
            # Test that different airfoil_ids get different polar data
            @test wing.sections[1].aero_model == POLAR_VECTORS
            @test wing.sections[2].aero_model == POLAR_VECTORS
            
            # Verify geometric progression along wingspan
            @test wing.sections[1].LE_point[2] == 5.0   # First section at y=5
            @test wing.sections[4].LE_point[2] == 0.0   # Middle section at y=0
            @test wing.sections[7].LE_point[2] == -5.0  # Last section at y=-5
        end
    end
    
    # Cleanup after all tests
    cleanup_test_files()
end
