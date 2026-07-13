using Test
using VortexStepMethod.ObjAdapter
using VortexStepMethod
using VortexStepMethod.AirfoilAero: NeuralFoilSolver
using LinearAlgebra

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

    @testset "generate_section_polars writes per-slice files" begin
        secdir = mktempdir()
        generate_section_polars(obj_path, secdir;
            n_slices=3, Re=5e5, alpha_range=-4:2:4, delta_range=-2:2:2,
            solver=NeuralFoilSolver(model_size="medium"), verbose=false)
        @test isfile(joinpath(secdir, "1.dat"))
        @test isfile(joinpath(secdir, "1.csv"))
        # a delta_range writes long-format POLAR_MATRICES + a .dat per deflection
        @test occursin("delta", lowercase(readline(joinpath(secdir, "1.csv"))))
        @test isfile(joinpath(secdir, "1_d2.dat"))

        @test_throws ErrorException generate_section_polars("missing.obj", secdir;
            n_slices=3, Re=5e5)
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
