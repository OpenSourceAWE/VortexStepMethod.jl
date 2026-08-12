using Test
using VortexStepMethod.ObjAdapter
using VortexStepMethod
using VortexStepMethod.AirfoilAero: NeuralFoilSolver
using LinearAlgebra
import YAML

obj_path = normpath(joinpath(@__DIR__, "..", "..",
                             "data", "ram_air_kite", "ram_air_kite_body.obj"))

@testset "ObjAdapter" begin
    @assert isfile(obj_path) "test obj mesh missing: $obj_path"

    # Reuse the suite-wide generated matrix geometry (default config, keyed and
    # cached under test/generated/) so the slow NeuralFoil conversion is shared with
    # the other tests that call ram_air_matrix_wing rather than regenerated here.
    out, yaml = ram_air_matrix_dir()

    @testset "perpendicular_sections slices the mesh" begin
        vertices, faces = read_faces(obj_path)
        secs = perpendicular_sections(vertices, faces, 4)
        @test length(secs) == 4
        for s in secs
            @test length(s.LE_point) == 3
            @test length(s.TE_point) == 3
            @test norm(s.TE_point .- s.LE_point) > 0
            @test !isempty(s.x_airfoil)
        end
    end

    @testset "obj_to_yaml (alpha,delta) matrices -> loadable Wing (NeuralFoil)" begin
        @test isfile(yaml)
        @test isfile(joinpath(out, "polars", "1.csv"))
        @test isfile(joinpath(out, "airfoils", "1.dat"))
        # a delta column marks the CSV as long-format POLAR_MATRICES
        @test occursin("delta", lowercase(readline(joinpath(out, "polars", "1.csv"))))

        wing = Wing(yaml; n_panels=6)
        body_aero = BodyAerodynamics([wing])
        @test length(body_aero.panels) == 6
    end

    @testset "generated geometry exposes its airfoils" begin
        airfoils = airfoils_from_yaml(yaml)
        @test !isempty(airfoils)
        @test all(af -> !isempty(af.x), airfoils)
    end

    @testset "obj_to_yaml writes shape + polar + Cp/cf per airfoil" begin
        outdir = mktempdir()
        yaml = obj_to_yaml(obj_path, outdir; n_sections=3, Re=5e5,
            alpha_range=-4:2:4, aero_solver=NeuralFoilSolver(model_size="medium"),
            verbose=false)
        @test isfile(yaml)
        info = Dict(YAML.load_file(yaml)["wing_airfoils"]["data"][1][3])
        @test all(haskey(info, k) for k in ("dat_file", "csv_file_path", "cp_file", "cf_file"))
        @test isfile(joinpath(outdir, info["dat_file"]))
        @test isfile(joinpath(outdir, info["cp_file"]))
        @test isfile(joinpath(outdir, info["cf_file"]))

        @test_throws ErrorException obj_to_yaml("missing.obj", outdir; n_sections=3, Re=5e5)
    end

    @testset "obj_to_yaml migrates existing node tables to another format" begin
        outdir = mktempdir()
        csv_yaml = obj_to_yaml(obj_path, outdir; n_sections=3, Re=5e5,
            alpha_range=-4:2:4, aero_solver=NeuralFoilSolver(model_size="medium"),
            verbose=false)
        csv_info = Dict(YAML.load_file(csv_yaml)["wing_airfoils"]["data"][1][3])

        arrow_yaml = obj_to_yaml(obj_path, outdir; n_sections=3, Re=5e5,
            verbose=false, table_format=:arrow)
        info = Dict(YAML.load_file(arrow_yaml)["wing_airfoils"]["data"][1][3])
        @test endswith(info["cp_file"], ".arrow")
        @test endswith(info["cf_file"], ".arrow")
        @test isfile(joinpath(outdir, info["cp_file"]))
        @test isfile(joinpath(outdir, info["cf_file"]))
        @test isfile(joinpath(outdir, csv_info["cp_file"]))
        @test Wing(arrow_yaml; n_panels=4) isa Wing

        # already in that format: nothing to convert, YAML untouched
        again = obj_to_yaml(obj_path, outdir; n_sections=3, Re=5e5,
            verbose=false, table_format=:arrow)
        @test Dict(YAML.load_file(again)["wing_airfoils"]["data"][1][3]) == info

        @test_throws ArgumentError ObjAdapter.migrate_node_tables(arrow_yaml, outdir,
                                                                 :parquet)
        rm(joinpath(outdir, info["cp_file"]))
        @test_throws ErrorException ObjAdapter.migrate_node_tables(arrow_yaml, outdir,
                                                                  :csv)
    end

    @testset "center_to_com! rejects non-triangular faces" begin
        verts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
        @test_throws ArgumentError center_to_com!(verts, [[1, 2, 3, 4]]; prn=false)
    end

    @testset "write_yaml emits nested and scalar values" begin
        dir = mktempdir()
        nested = joinpath(dir, "nested.yaml")
        write_yaml(nested, Dict("a" => [Dict("b" => Dict("c" => [1.0, 2.0]))]))
        text = read(nested, String)
        @test occursin("a:", text)
        @test occursin("c:", text)

        scalar = joinpath(dir, "scalar.yaml")
        write_yaml(scalar, 42.0)
        @test occursin("42", read(scalar, String))
    end
end
