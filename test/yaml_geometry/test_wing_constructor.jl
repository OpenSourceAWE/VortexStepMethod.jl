using VortexStepMethod
using VortexStepMethod: load_polar_data
using LinearAlgebra
using Test
using YAML
using Logging

@testset "Wing Constructor Tests" begin
    # Setup temporary files for testing
    test_yaml_path = joinpath(tempdir(), "test_wing.yaml")
    test_polar_dir = joinpath(tempdir(), "polars")
    
    # Clean up function
    function cleanup_test_files()
        test_dir = dirname(test_yaml_path)
        for file in [test_yaml_path, 
                     joinpath(test_dir, "standard_airfoil.csv"),
                     joinpath(test_dir, "alternate_airfoil.csv")]
            isfile(file) && rm(file; force=true)
        end
        isdir(test_polar_dir) && rm(test_polar_dir; recursive=true, force=true)
    end
    
    # Create polar data directory and files
    mkpath(test_polar_dir)
    
    # Copy the actual polar files to the temp directory for tests that reference them
    cp(test_data_path("yaml_geometry", "standard_airfoil.csv"), joinpath(test_polar_dir, "1.csv"); force=true)
    cp(test_data_path("yaml_geometry", "alternate_airfoil.csv"), joinpath(test_polar_dir, "2.csv"); force=true)
    
    # Also copy them to the test directory itself for direct reference tests
    test_dir = dirname(test_yaml_path)
    cp(test_data_path("yaml_geometry", "standard_airfoil.csv"), joinpath(test_dir, "standard_airfoil.csv"); force=true)
    cp(test_data_path("yaml_geometry", "alternate_airfoil.csv"), joinpath(test_dir, "alternate_airfoil.csv"); force=true)
    
    @testset "Valid YAML Wing Construction" begin
        # Use the actual YAML file from the test data
        cp(test_data_path("yaml_geometry", "simple_wing.yaml"), test_yaml_path; force=true)
        
        wing = Wing(test_yaml_path; n_panels=4)

        @test wing isa Wing
        @test wing.n_panels == 4
        @test wing.n_unrefined_sections == 2
        @test wing.spanwise_distribution == LINEAR
        @test wing.spanwise_direction ≈ [0.0, 1.0, 0.0]
        @test length(wing.unrefined_sections) == 2  # simple_wing has 2 sections
        
        # Test section coordinates (sections are sorted by spanwise position)
        # simple_wing.yaml has: [1, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0] and [1, 0.0, -1.0, 0.0, 1.0, -1.0, 0.0]
        # Sorted by y-coordinate: y=1.0, y=-1.0
        @test wing.unrefined_sections[1].LE_point ≈ [0.0, 1.0, 0.0]
        @test wing.unrefined_sections[1].TE_point ≈ [1.0, 1.0, 0.0]
        @test wing.unrefined_sections[2].LE_point ≈ [0.0, -1.0, 0.0]
        @test wing.unrefined_sections[2].TE_point ≈ [1.0, -1.0, 0.0]
        
        # Test that sections have polar data
        @test wing.unrefined_sections[1].aero_model == POLAR_VECTORS
        @test wing.unrefined_sections[2].aero_model == POLAR_VECTORS
        
        # Test polar data is loaded
        @test wing.unrefined_sections[1].aero_data isa Tuple
        @test length(wing.unrefined_sections[1].aero_data) == 4
        @test length(wing.unrefined_sections[1].aero_data[1]) >= 3  # at least 3 alpha points
    end
    
    @testset "Wing Constructor Parameters" begin
        # Use the simple_wing.yaml and test custom parameters
        cp(test_data_path("yaml_geometry", "simple_wing.yaml"), test_yaml_path; force=true)
        
        # Test custom parameters
        wing = Wing(
            test_yaml_path;
            n_panels=8,
            spanwise_distribution=COSINE,
            remove_nan=false
        )

        @test wing.n_panels == 8
        @test wing.n_unrefined_sections == 2
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
    - [1, 0.0, -1.0, 0.0, 1.0, -1.0, 0.0]

wing_airfoils:
  alpha_range: [-10, 10, 1]
  reynolds: 1000000
  headers: [airfoil_id, type, info_dict]
  data:
    - [1, polars, {csv_file_path: "nonexistent.csv"}]
"""
        write(test_yaml_path, yaml_content)
        
        wing = suppress_warnings(() -> Wing(test_yaml_path; n_panels=2))
        
        # Should fall back to INVISCID model
        @test wing.unrefined_sections[1].aero_model == INVISCID
        @test wing.unrefined_sections[1].aero_data === nothing
    end
    
    @testset "Sections Without Polar Files" begin
        # Create YAML with empty polar file paths
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
    - [1, polars, {csv_file_path: ""}]
"""
        write(test_yaml_path, yaml_content)
        
        wing = suppress_warnings(() -> Wing(test_yaml_path; n_panels=2))
        
        # Should fall back to INVISCID model
        @test wing.unrefined_sections[1].aero_model == INVISCID
    end
    
    @testset "Invalid Parameters" begin
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

        # Test # n_groups=0 (no grouping functionality) - backward compatibility
        wing_no_groups = Wing(test_yaml_path; n_panels=4)
        @test wing_no_groups.n_panels == 4

        # Test invalid spanwise direction
        @test_throws ArgumentError Wing(test_yaml_path; spanwise_direction=[1.0, 0.0, 0.0])
    end
    
    @testset "Relative Path Resolution" begin
        # Test that relative paths in CSV files are resolved relative to YAML file
        subdir = joinpath(tempdir(), "subtest")
        mkpath(subdir)
        
        # Copy the simple wing file to subdirectory
        subdir_yaml_path = joinpath(subdir, "wing.yaml")
        cp(test_data_path("yaml_geometry", "simple_wing.yaml"), subdir_yaml_path; force=true)
        
        # Copy the polar file to subdirectory too
        subdir_polar_path = joinpath(subdir, "standard_airfoil.csv")
        cp(test_data_path("yaml_geometry", "standard_airfoil.csv"), subdir_polar_path; force=true)
        
        wing = Wing(subdir_yaml_path; n_panels=2)
        
        # Should successfully load polar data with relative path
        @test wing.unrefined_sections[1].aero_model == POLAR_VECTORS
        @test wing.unrefined_sections[1].aero_data isa Tuple
        @test wing.unrefined_sections[2].aero_model == POLAR_VECTORS
        @test wing.unrefined_sections[2].aero_data isa Tuple
        
        # Cleanup
        rm(subdir; recursive=true)
    end
    
    @testset "Complex Wing Geometry" begin
        # Use the actual complex_wing.yaml file
        cp(test_data_path("yaml_geometry", "complex_wing.yaml"), test_yaml_path; force=true)

        wing = Wing(test_yaml_path; n_panels=12)

        @test wing.n_panels == 12
        @test wing.n_unrefined_sections == 7
        @test length(wing.unrefined_sections) == 7
        
        # Test that different airfoil_ids get different polar data
        @test wing.unrefined_sections[1].aero_model == POLAR_VECTORS
        @test wing.unrefined_sections[2].aero_model == POLAR_VECTORS
        
        # Verify geometric progression along wingspan
        @test wing.unrefined_sections[1].LE_point[2] == 5.0   # First section at y=5
        @test wing.unrefined_sections[4].LE_point[2] == 0.0   # Middle section at y=0
        @test wing.unrefined_sections[7].LE_point[2] == -5.0  # Last section at y=-5
    end
    
    @testset "VSMSettings Constructor" begin
        # Test Wing constructor that takes VSMSettings
        # Use the actual simple_wing.yaml file
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        
        # Create VSMSettings with wing configuration
        settings = VSMSettings()
        settings.wings = [WingSettings(
            geometry_file=simple_wing_file,
            n_panels=6,
            spanwise_panel_distribution=COSINE
        )]

        # Test Wing constructor with VSMSettings
        wing = Wing(settings)

        @test wing isa Wing
        @test wing.n_panels == 6
        @test wing.n_unrefined_sections == 2
        @test wing.spanwise_distribution == COSINE
        @test length(wing.unrefined_sections) == 2
        @test wing.unrefined_sections[1].aero_model == POLAR_VECTORS
        @test wing.unrefined_sections[2].aero_model == POLAR_VECTORS
    end
    
    @testset "Shared Test Data Usage" begin
        # Demonstrate using module-specific test data files
        simple_wing_file = test_data_path("yaml_geometry", "simple_wing.yaml")
        @test isfile(simple_wing_file)
        
        # Test basic Wing construction with shared data
        wing = Wing(simple_wing_file; n_panels=4)
        @test wing isa Wing
        @test wing.n_panels == 4
        @test wing.n_unrefined_sections == 2
        @test length(wing.unrefined_sections) == 2
        
        # Test complex wing construction
        complex_wing_file = test_data_path("yaml_geometry", "complex_wing.yaml")
        @test isfile(complex_wing_file)
        
        complex_wing = Wing(complex_wing_file; n_panels=12)
        @test complex_wing isa Wing
        @test complex_wing.n_panels == 12
        @test complex_wing.n_unrefined_sections == 7
        @test length(complex_wing.unrefined_sections) == 7
        
        # Verify polar data is loaded from shared files
        @test complex_wing.unrefined_sections[1].aero_model == POLAR_VECTORS
        @test complex_wing.unrefined_sections[1].aero_data isa Tuple
        
        # Test with module-specific convenience function - create a standard wing for this test
        standard_wing_file = simple_wing_file  # Use simple_wing as our "standard"
        @test isfile(standard_wing_file)
        
        standard_wing = Wing(standard_wing_file; n_panels=2)
        @test standard_wing isa Wing
        @test length(standard_wing.unrefined_sections) == 2
    end
    
    # Cleanup after all tests
    cleanup_test_files()
end
