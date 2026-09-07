using VortexStepMethod
using VortexStepMethod: ObjAdapter
using Test

ram_air_dir = joinpath(dirname(dirname(@__DIR__)), "data", "ram_air_kite")

@testset "Test settings.jl" begin
    vss = VSMSettings(joinpath(ram_air_dir, "vsm_settings_dual.yaml"))
    @test vss isa VSMSettings
    @test vss.solver_settings isa SolverSettings
    @test vss.wings isa Vector{WingSettings}
    @test length(vss.wings) == 2
    io = IOBuffer(repr(vss))
    @test countlines(io) > 0

    # A file that names neither block still loads, on the defaults.
    bare = vss.wings[1]
    @test bare.mesh isa MeshSettings
    @test bare.airfoil isa AirfoilSettings
    @test isnothing(delta_range(bare.airfoil))
    @test airfoil_solver(bare.airfoil) isa
        VortexStepMethod.AirfoilAero.NeuralFoilSolver
end

@testset "mesh and airfoil blocks" begin
    yaml = """
    wings:
      - name: wing
        n_panels: 4
        spanwise_panel_distribution: LINEAR
        spanwise_direction: [0.0, 1.0, 0.0]
        remove_nan: true
        mesh:
          obj_file: obj/wing.obj
          n_sections: 45
          rotation: [[0, 0, -1], [-1, 0, 0], [0, 1, 0]]
          clearance: 0.0
        airfoil:
          solver: xfoil
          n_crit: 4.0
          alpha_range: [-15, 3, 90]
          delta_range: [-40, 10, 40]
          v_app: 25.0
          chord_ref: 6.0
          table_format: csv
    solver_settings:
      aerodynamic_model_type: VSM
      type_initial_gamma_distribution: ELLIPTIC
      density: 1.225
      mu: 1.81e-5
    """
    path = tempname() * ".yaml"
    write(path, yaml)
    set = VSMSettings(path; data_prefix = false)
    wing = set.wings[1]

    @test wing.mesh.n_sections == 45
    # the mesh the sections come from, not the wing's own obj_file, which is an
    # alternative to geometry_file and cannot be given beside one
    @test wing.mesh.obj_file == "obj/wing.obj"
    @test isempty(wing.obj_file)
    @test rotation_matrix(wing.mesh) == [0 0 -1; -1 0 0; 0 1 0]
    @test alpha_range(wing.airfoil) == -15.0:3.0:90.0
    @test delta_range(wing.airfoil) == -40.0:10.0:40.0
    # the air is stated once, in solver_settings, and Reynolds reads it there
    @test reynolds(set, wing) ≈ 1.225 * 25.0 * 6.0 / 1.81e-5
    @test airfoil_solver(wing.airfoil) isa VortexStepMethod.AirfoilAero.XFoilSolver
    # tables and live polars come off the same block, so they cannot disagree
    live = VortexStepMethod.AirfoilAero.LivePolarSettings(wing.airfoil)
    @test live.n_crit == wing.airfoil.n_crit
    @test live.model_size == wing.airfoil.model_size
    @test keys(slice_args(wing)) ==
        (:rotation, :n_bins, :wrap_method, :wingtip_distance)
    @test preview_args(wing).n_slices == 45
    @test preview_args(wing).n_bins == wing.mesh.n_bins
    # the generator takes a Symbol, so the settings carry one
    @test wing.airfoil.table_format === :csv

    @test_throws ArgumentError VortexStepMethod.airfoil_settings(
        Dict("solver" => "rans"))
    @test_throws ArgumentError VortexStepMethod.airfoil_settings(
        Dict("table_format" => "parquet"))
    @test_throws Exception VortexStepMethod.airfoil_settings(Dict("n_crt" => 4.0))
end

@testset "an unnamed mesh block slices as an unconfigured call" begin
    vertices, faces = ObjAdapter.read_faces(joinpath(ram_air_dir,
                                                    "ram_air_kite_body.obj"))
    args = slice_args(WingSettings())

    @test ObjAdapter.perpendicular_sections(vertices, faces, 4; args.rotation,
              args.n_bins, args.wingtip_distance) ==
          ObjAdapter.perpendicular_sections(vertices, faces, 4)
    @test args.wrap_method == VortexStepMethod.AirfoilAero.ShrinkWrap()
end
nothing
