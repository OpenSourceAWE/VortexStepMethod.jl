using VortexStepMethod
using VortexStepMethod: calculate_cl, calculate_cd_cm, calculate_projected_area, calculate_AIC_matrices!, reinit!
using DifferentiationInterface
using LinearAlgebra
using Test
using Logging

# Helper to find the data directory for test files
function _find_ram_data_dir()
    data_dir = joinpath(dirname(dirname(@__DIR__)), "data", "ram_air_kite")
    if isdir(data_dir)
        return data_dir
    else
        # Fallback for case where test is run from different working directory
        return joinpath(@__DIR__, "..", "..", "data", "ram_air_kite")
    end
end

if !@isdefined ram_wing_results
    data_dir = _find_ram_data_dir()
    body_path = joinpath(tempdir(), "ram_air_kite_body.obj")
    foil_path = joinpath(tempdir(), "ram_air_kite_foil.dat")

    body_src = joinpath(data_dir, "ram_air_kite_body.obj")
    foil_src = joinpath(data_dir, "ram_air_kite_foil.dat")

    if isfile(body_src) && isfile(foil_src)
        cp(body_src, body_path; force=true)
        cp(foil_src, foil_path; force=true)
        for kind in ("cl", "cd", "cm")
            name = "ram_air_kite_foil_$(kind)_polar.csv"
            cp(joinpath(data_dir, name), joinpath(tempdir(), name); force=true)
        end
    else
        error("Required data files not found: $body_src or $foil_src")
    end

    ram_wing = ram_air_matrix_wing(; n_panels=8, n_sections=4,
        alpha_range=deg2rad.(-5:5:15),
        delta_range=deg2rad.(-3:3:3),
    )
end

@testset "Linearize Jacobian validation" begin
    # The Jacobian is computed by FD with absstep=relstep=fd_step inside
    # linearize. Two assertions per perturbation direction:
    #   1. perturb by exactly the FD step → predicted change must match the
    #      observed change to (near) numerical precision; this verifies the
    #      Jacobian is internally consistent with its own FD computation.
    #   2. perturb by 2× the FD step → must still match within a looser
    #      tolerance dominated by the O(h²) Taylor remainder; this verifies
    #      the Jacobian columns reflect true local sensitivity rather than
    #      numerical noise (a noise-driven Jacobian would not extrapolate).

    va = [15.0, 1.0, 0.5]
    theta = deg2rad.([2.0, 1.0, -1.0, -2.0])
    delta = deg2rad.([1.0, 0.5, -0.5, -1.0])
    omega = [0.0, 0.1, 0.0]

    fd_step = 1e-3

    VortexStepMethod.unrefined_deform!(ram_wing, theta, delta; smooth=false)
    body_aero = BodyAerodynamics([ram_wing]; va, omega)
    solver = Solver(body_aero;
        aerodynamic_model_type=VSM,
        is_with_artificial_damping=false,
        atol=1e-5,
        rtol=1e-5,
        solver_type=NONLIN,
    )

    base_inputs = [theta; va; omega; delta]
    jac, lin_res, lin_converged = VortexStepMethod.linearize(
        solver, body_aero, base_inputs;
        theta_idxs=1:4, va_idxs=5:7, omega_idxs=8:10, delta_idxs=11:14,
        moment_frac=0.1,
        backend=AutoFiniteDiff(absstep=fd_step, relstep=fd_step),
    )
    @test lin_converged

    # Warm-start to match the protocol linearize uses internally; with a
    # finite-difference Jacobian Newton, cold-start vs warm-start can
    # differ at the noise floor (~sqrt(eps)) and that floor dominates
    # the small Δ used at scale=1.
    function evaluate_at!(input_vec)
        perturbed_theta = input_vec[1:4]
        perturbed_va    = input_vec[5:7]
        perturbed_omega = input_vec[8:10]
        perturbed_delta = input_vec[11:14]
        VortexStepMethod.unrefined_deform!(
            ram_wing, perturbed_theta, perturbed_delta; smooth=false)
        reinit!(body_aero; init_aero=false,
            va=perturbed_va, omega=perturbed_omega)
        VortexStepMethod.solve!(solver, body_aero; log=false)
        return [solver.sol.force; solver.sol.moment;
                solver.sol.moment_unrefined_dist]
    end

    baseline_res = evaluate_at!(base_inputs)
    @test baseline_res ≈ lin_res rtol=1e-5

    @testset "scale=$scale" for scale in (1.0, 2.0)
        for j in 1:length(base_inputs)
            step = max(fd_step, fd_step * abs(base_inputs[j])) * scale
            perturbed = copy(base_inputs)
            perturbed[j] += step

            actual_change = evaluate_at!(perturbed) .- baseline_res
            predicted_change = jac[:, j] .* step

            scale_norm = max(norm(actual_change), 1e-8)
            err_rel = norm(predicted_change .- actual_change) / scale_norm
            tol = scale == 1.0 ? 1e-3 : 1e-1
            @test err_rel < tol
        end
    end
end
