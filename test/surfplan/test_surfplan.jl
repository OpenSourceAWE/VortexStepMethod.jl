using Test
using VortexStepMethod
using VortexStepMethod.SurfplanAdapter: surfplan_to_aero_yaml
using VortexStepMethod.AirfoilAero: NeuralFoilSolver
import YAML

@testset "SurfplanAdapter" begin
    dat_src = joinpath(@__DIR__, "..", "airfoil_aero", "data", "test_airfoil.dat")
    @assert isfile(dat_src) "missing test airfoil: $dat_src"

    adapter_dir = mktempdir()
    output_dir = mktempdir()
    cp(dat_src, joinpath(adapter_dir, "airfoil.dat"))

    aero = Dict(
        "wing_airfoils" => Dict(
            "reynolds" => 5.0e5,
            "headers" => ["id", "name", "info_dict"],
            "data" => Any[Any[1, "test", Dict("dat_file_path" => "airfoil.dat")]],
        ),
        "wing_sections" => Dict(
            "headers" => ["airfoil_id", "LE_x", "LE_y", "LE_z",
                          "TE_x", "TE_y", "TE_z"],
            "data" => Any[Any[1, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                          Any[1, 0.0, 2.0, 0.0, 1.0, 2.0, 0.0]],
        ),
    )
    YAML.write_file(joinpath(adapter_dir, "aero_geometry.yaml"), aero)

    yaml_path = surfplan_to_aero_yaml(adapter_dir, output_dir;
        aero_solver=NeuralFoilSolver(model_size="medium"),
        alpha_range=-4:4:4, Re=5.0e5, verbose=false)

    @test isfile(yaml_path)
    @test isfile(joinpath(output_dir, "airfoils", "1.dat"))
    @test isfile(joinpath(output_dir, "polars", "1.csv"))

    wing = Wing(yaml_path)
    @test wing isa Wing

    # a missing dat_file_path is a hard error (the profiles are required)
    aero["wing_airfoils"]["data"][1][3] = Dict("dat_file_path" => "")
    YAML.write_file(joinpath(adapter_dir, "aero_geometry.yaml"), aero)
    @test_throws ErrorException surfplan_to_aero_yaml(adapter_dir, mktempdir();
        alpha_range=-4:4:4, Re=5.0e5, verbose=false)
end
