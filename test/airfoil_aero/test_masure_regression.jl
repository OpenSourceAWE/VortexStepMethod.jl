using Test
import YAML
using VortexStepMethod
using VortexStepMethod.ObjAdapter: resolve_aero_geometry
using VortexStepMethod.AirfoilAero: MASURE_PARAMETERS, load_masure_model, masure_aero,
                                    write_yaml, npzread

fixture_dir = joinpath(@__DIR__, "data", "masure")

@testset "masure regression" begin
    @testset "Extra-Trees evaluation matches sklearn predict, on-split rows too" begin
        reference = npzread(joinpath(fixture_dir, "reference.npz"))
        model = load_masure_model(1e6, fixture_dir)
        for (row, expected) in zip(eachrow(reference["X"]), eachrow(reference["Y"]))
            params = Dict(zip(MASURE_PARAMETERS, row[1:6]))
            cl, cd, cm = masure_aero(model, params, [row[7]])
            @test [cd[1], cl[1], cm[1]] ≈ expected atol=1e-12
        end
    end

    @testset "load_masure_model rejects an untrained Re and a missing file" begin
        @test_throws "No masure regression model for Re = 3.0e6" load_masure_model(
            3e6, fixture_dir)
        @test_throws "Masure regression model not found" load_masure_model(
            5e6, fixture_dir)
    end

    @testset "resolve_aero_geometry turns masure_regression into a loadable polar" begin
        dir = mktempdir()
        yaml_in = joinpath(dir, "awesio.yaml")
        params = Dict("t" => 0.08, "eta" => 0.2, "kappa" => 0.09, "delta" => -2.0,
                      "lambda" => 0.2, "phi" => 0.6)
        write_yaml(yaml_in, Dict(
            "wing_sections" => Dict(
                "headers" => ["airfoil_id", "LE_x", "LE_y", "LE_z", "TE_x", "TE_y", "TE_z"],
                "data" => [Any[1, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0],
                           Any[1, 0.0, -1.0, 0.0, 1.0, -1.0, 0.0]]),
            "wing_airfoils" => Dict(
                "alpha_range" => [-4, 4, 2], "reynolds" => 1e6,
                "headers" => ["airfoil_id", "type", "info_dict"],
                "data" => [Any[1, "masure_regression", params]])))
        @test_throws "needs ml_models_dir" resolve_aero_geometry(
            yaml_in, joinpath(dir, "out"); verbose=false)
        yaml_out = resolve_aero_geometry(yaml_in, joinpath(dir, "out");
                                         ml_models_dir=fixture_dir, verbose=false)
        airfoil = YAML.load_file(yaml_out)["wing_airfoils"]["data"][1]
        @test airfoil[2] == "polars"
        section = Wing(yaml_out; n_panels=2).unrefined_sections[1]
        @test section.aero_model == POLAR_VECTORS
        cl, _, _ = masure_aero(load_masure_model(1e6, fixture_dir), params, [2.0])
        @test section.aero_data[1][4] ≈ deg2rad(2.0)
        @test section.aero_data[2][4] ≈ cl[1]
    end
end
