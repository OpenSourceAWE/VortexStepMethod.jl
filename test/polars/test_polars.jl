using VortexStepMethod
using AirfoilAero
using Test
using DelimitedFiles

@testset "XFoil polar generation" begin
    data_dir = joinpath(dirname(dirname(@__DIR__)), "data", "ram_air_kite")
    foil_src = joinpath(data_dir, "ram_air_kite_foil.dat")
    @assert isfile(foil_src) "fixture foil file missing: $foil_src"

    work_dir = mktempdir()
    foil_path = joinpath(work_dir, "test_foil.dat")
    cp(foil_src, foil_path; force=true)

    cl_path = joinpath(work_dir, "test_foil_cl_polar.csv")
    cd_path = joinpath(work_dir, "test_foil_cd_polar.csv")
    cm_path = joinpath(work_dir, "test_foil_cm_polar.csv")

    alpha_range = deg2rad.(-2:1:2)
    delta_range = deg2rad.(-1:1:1)

    AirfoilAero.create_polars(;
        dat_path=foil_path,
        cl_polar_path=cl_path,
        cd_polar_path=cd_path,
        cm_polar_path=cm_path,
        wind_vel=15.0,
        area=20.0,
        width=8.0,
        crease_frac=0.75,
        alpha_range,
        delta_range,
    )

    @test isfile(cl_path)
    @test isfile(cd_path)
    @test isfile(cm_path)

    cl_matrix, cl_alpha, cl_delta = VortexStepMethod.read_aero_matrix(cl_path)
    @test size(cl_matrix) == (length(alpha_range), length(delta_range))
    @test cl_alpha ≈ collect(alpha_range)
    @test cl_delta ≈ collect(delta_range)
    @test all(isfinite, filter(!isnan, cl_matrix))
    @test count(isnan, cl_matrix) <= length(cl_matrix) ÷ 4

    cd_matrix, _, _ = VortexStepMethod.read_aero_matrix(cd_path)
    @test size(cd_matrix) == (length(alpha_range), length(delta_range))
    @test all(x -> isnan(x) || x > 0, cd_matrix)

    cm_matrix, _, _ = VortexStepMethod.read_aero_matrix(cm_path)
    @test size(cm_matrix) == (length(alpha_range), length(delta_range))
end
