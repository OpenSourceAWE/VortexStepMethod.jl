using Test
using VortexStepMethod.ObjAdapter
using VortexStepMethod
using VortexStepMethod: load_polar_data
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
        @test all(haskey(info, k)
                  for k in ("dat_file", "polar_file_path", "cp_file", "cf_file"))
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
        write(csv_yaml, replace(read(csv_yaml, String),
                                "polar_file_path" => "csv_file_path"))

        arrow_yaml = obj_to_yaml(obj_path, outdir; n_sections=3, Re=5e5,
            verbose=false, table_format=:arrow)
        info = Dict(YAML.load_file(arrow_yaml)["wing_airfoils"]["data"][1][3])
        @test !haskey(info, "csv_file_path")
        for key in ("cp_file", "cf_file", "polar_file_path")
            @test endswith(info[key], ".arrow")
            @test isfile(joinpath(outdir, info[key]))
            @test isfile(joinpath(outdir, csv_info[key]))
        end
        @test isequal(load_polar_data(joinpath(outdir, info["polar_file_path"])),
                      load_polar_data(joinpath(outdir, csv_info["polar_file_path"])))
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

    @testset "geometry_path writes the YAML outside the table directory" begin
        root = mktempdir()
        tables = joinpath(root, "tables")
        yaml_path = joinpath(root, "nf_aero_geometry.yaml")
        written = obj_to_yaml(obj_path, tables; n_sections=3, Re=5e5,
            verbose=false, geometry_path=yaml_path)
        @test written == yaml_path
        @test isfile(yaml_path)
        @test !isfile(joinpath(tables, "geometry.yaml"))

        info = Dict(YAML.load_file(yaml_path)["wing_airfoils"]["data"][1][3])
        for key in ("polar_file_path", "dat_file")
            @test startswith(info[key], "tables/")
            @test isfile(joinpath(dirname(yaml_path), info[key]))
        end
        # The loader resolves references against the YAML's own directory, so the
        # prefixed paths have to be what makes this work.
        @test Wing(yaml_path; n_panels=4) isa Wing

        @test ObjAdapter.table_path_prefix(joinpath(tables, "geometry.yaml"),
                                           tables) == ""
    end

    @testset "resolve_aero_geometry writes neuralfoil polars as table_format" begin
        dir = mktempdir()
        cp(joinpath(@__DIR__, "..", "airfoil_aero", "data", "test_airfoil.dat"),
           joinpath(dir, "airfoil.dat"))
        yaml_in = joinpath(dir, "awesio.yaml")
        write_yaml(yaml_in, Dict(
            "wing_sections" => Dict(
                "headers" => ["airfoil_id", "LE_x", "LE_y", "LE_z", "TE_x", "TE_y", "TE_z"],
                "data" => [Any[1, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0],
                           Any[1, 0.0, -1.0, 0.0, 1.0, -1.0, 0.0]]),
            "wing_airfoils" => Dict(
                "alpha_range" => [-4, 4, 2], "reynolds" => 5e5,
                "headers" => ["airfoil_id", "type", "info_dict"],
                "data" => [Any[1, "neuralfoil", Dict("dat_file_path" => "airfoil.dat",
                                                     "model_size" => "medium")]])))
        yaml_out = resolve_aero_geometry(yaml_in, joinpath(dir, "out");
                                         table_format=:arrow, verbose=false)
        info = YAML.load_file(yaml_out)["wing_airfoils"]["data"][1][3]
        @test endswith(info["polar_file_path"], "1.arrow")
        @test isfile(info["polar_file_path"])
        @test Wing(yaml_out; n_panels=2).unrefined_sections[1].aero_model == POLAR_VECTORS
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
