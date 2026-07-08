using Test
using VortexStepMethod.ObjAdapter
using VortexStepMethod
using VortexStepMethod.AirfoilAero: NeuralFoilSolver
using LinearAlgebra

obj_path = normpath(joinpath(@__DIR__, "..", "..",
                             "data", "ram_air_kite", "ram_air_kite_body.obj"))

@testset "ObjAdapter" begin
    @assert isfile(obj_path) "test obj mesh missing: $obj_path"

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
        out = mktempdir()
        yaml = obj_to_yaml(obj_path, out;
            n_sections=3, Re=5e5,
            alpha_range=-4:2:4, delta_range=-2:2:2,
            aero_solver=NeuralFoilSolver(), verbose=false)
        @test isfile(yaml)
        @test isfile(joinpath(out, "polars", "1.csv"))
        @test isfile(joinpath(out, "airfoils", "1.dat"))
        # each trailing-edge deflection is stored as its own .dat
        @test isfile(joinpath(out, "airfoils", "1_d2.dat"))
        @test isfile(joinpath(out, "airfoils", "1_d-2.dat"))
        # a delta column marks the CSV as long-format POLAR_MATRICES
        @test occursin("delta", lowercase(readline(joinpath(out, "polars", "1.csv"))))

        wing = Wing(yaml; n_panels=6)
        body_aero = BodyAerodynamics([wing])
        @test length(body_aero.panels) == 6
    end

    @testset "obj_to_yaml -> per-section polars (NeuralFoil)" begin
        out = mktempdir()
        yaml = obj_to_yaml(obj_path, out;
            n_sections=3, Re=5e5, alpha_range=-6:2:6,
            aero_solver=NeuralFoilSolver(model_size="large"), verbose=false)
        @test isfile(yaml)
        airfoils = airfoils_from_yaml(yaml)
        @test !isempty(airfoils)
        @test all(af -> !isempty(af.x), airfoils)

        wing = Wing(yaml; n_panels=6)
        @test length(BodyAerodynamics([wing]).panels) == 6
    end
end
