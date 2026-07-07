using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using GLMakie
using DifferentiationInterface
using LinearAlgebra
using VortexStepMethod
using ObjAdapter
using VortexStepMethod: linearize, unrefined_deform!, reinit!

# Sweep each linearize input around the operating point and overlay the
# FiniteDiff and ForwardDiff tangents on the sweep curve.

n_unrefined = 4

# Convert-then-load: obj mesh -> shared NeuralFoil (alpha, delta) POLAR_MATRICES.
yaml = obj_to_matrix_yaml(
    joinpath("data", "ram_air_kite", "ram_air_kite_body.obj"), mktempdir();
    n_sections=n_unrefined, wind_vel=15.0,
    alpha_range=deg2rad.(-10:2:20), delta_range=deg2rad.(-5:5:5), verbose=false)
wing = Wing(yaml; n_panels=16)
body_aero = BodyAerodynamics([wing])

solver = Solver(body_aero;
    aerodynamic_model_type=VSM,
    is_with_artificial_damping=false,
    rtol=1e-7,
    solver_type=LOOP,
    use_gamma_prev=false,
)

v_a       = 15.0
aoa_deg   = 10.0
aoa_rad   = deg2rad(aoa_deg)
side_slip = 0.0
va_b_0    = [
    cos(aoa_rad) * cos(side_slip),
    sin(side_slip),
    sin(aoa_rad),
] * v_a
omega_b_0 = zeros(3)
theta_0   = zeros(n_unrefined)

theta_idxs = 1:n_unrefined
va_idxs    = (n_unrefined + 1):(n_unrefined + 3)
omega_idxs = (n_unrefined + 4):(n_unrefined + 6)
y0         = [theta_0; va_b_0; omega_b_0]

@info "Computing FiniteDiff Jacobian …"
t_fd = @elapsed begin
    jac_fd, x0_fd, conv_fd = linearize(
        solver, body_aero, y0;
        theta_idxs, va_idxs, omega_idxs,
        aero_coeffs=true,
        backend=AutoFiniteDiff(absstep=1e-5, relstep=1e-5),
    )
end
conv_fd || @warn "FiniteDiff linearize did not converge at operating point"

@info "Computing ForwardDiff Jacobian …"
t_fwd = @elapsed begin
    jac_fwd, x0_fwd, conv_fwd = linearize(
        solver, body_aero, y0;
        theta_idxs, va_idxs, omega_idxs,
        aero_coeffs=true,
        backend=AutoForwardDiff(),
    )
end
conv_fwd || @warn "ForwardDiff linearize did not converge at operating point"

x0 = x0_fd

@info "Operating point" CFx=x0[1] CFy=x0[2] CFz=x0[3] CM=x0[4:6]
@info "Timing" t_finitediff_s=t_fd t_forwarddiff_s=t_fwd speedup=t_fd / t_fwd

denom = max(maximum(abs, jac_fwd), eps())
rel_err_fwd_fd = maximum(abs.(jac_fwd .- jac_fd)) / denom
@info "AutoForwardDiff vs AutoFiniteDiff" rel_err=rel_err_fwd_fd

input_labels = [
    ["θ_$i" for i in 1:n_unrefined]...,
    "va_x", "va_y", "va_z",
    "ω_x",  "ω_y",  "ω_z",
]
output_labels = [
    "CFx", "CFy", "CFz",
    "CMx", "CMy", "CMz",
    ["cm_$i" for i in 1:n_unrefined]...,
]
n_inputs  = length(input_labels)
n_outputs = length(output_labels)

@assert size(jac_fd) == (n_outputs, n_inputs)
@assert size(jac_fwd) == (n_outputs, n_inputs)

input_scales = [
    fill(0.05, n_unrefined)...,    # θ [rad] : ±0.05 rad ≈ ±2.9°
    fill(1.0, 3)...,               # va: ±1 m/s
    fill(0.05, 3)...,              # ω : ±0.05 rad/s
]
n_sweep      = 11
sweep_frac   = range(-1.0, 1.0; length=n_sweep)

results_sw = [zeros(n_sweep, n_outputs) for _ in 1:n_inputs]

last_theta = fill(NaN, n_unrefined)

function solve_at!(y)
    theta = y[theta_idxs]
    va    = y[va_idxs]
    omega = y[omega_idxs]
    if !all(theta .== last_theta)
        unrefined_deform!(wing, theta, nothing; smooth=false)
        reinit!(body_aero; init_aero=false)
        last_theta .= theta
    end
    set_va!(body_aero, va, omega)
    solve!(solver, body_aero; log=false)
    return [
        solver.sol.force_coeffs...,
        solver.sol.moment_coeffs...,
        solver.sol.cm_unrefined_dist...,
    ]
end

@info "Sweeping each input …"
for ci in 1:n_inputs
    @info "  $(input_labels[ci])"
    for (si, frac) in enumerate(sweep_frac)
        y = copy(y0)
        y[ci] += frac * input_scales[ci]
        results_sw[ci][si, :] .= solve_at!(y)
    end
end

solve_at!(y0)

fig = Figure(size=(180 * n_inputs + 80, 90 * n_outputs + 100))
Label(fig[0, 1:n_inputs],
    "linearize check: sweeps (blue) vs FiniteDiff tangent (red dashed) " *
    "vs ForwardDiff tangent (green dashed) — " *
    "FD: $(round(t_fd; digits=2)) s, FwdDiff: $(round(t_fwd; digits=2)) s, " *
    "speedup: $(round(t_fd / t_fwd; digits=2))×";
    fontsize=16, font=:bold, tellwidth=false)

for ri in 1:n_outputs
    for ci in 1:n_inputs
        ax = Axis(fig[ri, ci];
            xlabel = ri == n_outputs ? input_labels[ci] : "",
            ylabel = ci == 1         ? output_labels[ri] : "",
            xticklabelsvisible = ri == n_outputs,
            yticklabelsvisible = ci == 1,
            xticksvisible      = ri == n_outputs,
            yticksvisible      = ci == 1,
        )
        delta_input  = sweep_frac .* input_scales[ci]
        ys_sweep     = results_sw[ci][:, ri]
        ys_linear_fd  = x0[ri] .+ jac_fd[ri, ci]  .* delta_input
        ys_linear_fwd = x0[ri] .+ jac_fwd[ri, ci] .* delta_input

        lines!(ax, delta_input, ys_sweep;
            color=:steelblue, linewidth=1.5)
        lines!(ax, delta_input, ys_linear_fd;
            color=:crimson, linestyle=:dash, linewidth=1.2)
        lines!(ax, delta_input, ys_linear_fwd;
            color=:seagreen, linestyle=:dash, linewidth=1.2)
        scatter!(ax, [0.0], [x0[ri]];
            color=:black, markersize=6)
    end
end

colgap!(fig.layout, 6)
rowgap!(fig.layout, 4)

display(fig)

# Skip near-zero entries — relative error there is sweep rounding noise.
function worst_jac_vs_sweep(jac, label)
    sig_threshold = 1e-3 * maximum(abs, jac)
    max_rel = 0.0
    worst   = (0, 0, 0.0, 0.0)
    for ri in 1:n_outputs, ci in 1:n_inputs
        delta = sweep_frac .* input_scales[ci]
        mid   = div(n_sweep, 2) + 1
        slope = (results_sw[ci][mid + 1, ri] -
                 results_sw[ci][mid - 1, ri]) /
                (delta[mid + 1] - delta[mid - 1])
        jc = jac[ri, ci]
        max(abs(jc), abs(slope)) < sig_threshold && continue
        rel = abs(slope - jc) / max(abs(jc), abs(slope))
        if rel > max_rel
            max_rel = rel
            worst   = (ri, ci, slope, jc)
        end
    end
    @info "Worst $label vs sweep mismatch (significant entries)" rel=max_rel output=output_labels[worst[1]] input=input_labels[worst[2]] sweep_slope=worst[3] jac_entry=worst[4] sig_threshold
end

worst_jac_vs_sweep(jac_fd,  "FiniteDiff")
worst_jac_vs_sweep(jac_fwd, "ForwardDiff")

nothing
