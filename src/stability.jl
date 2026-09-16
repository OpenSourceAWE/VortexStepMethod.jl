"""
    stability_derivatives(solver, body_aero, alpha, beta, wind_speed; kwargs...)

Aerodynamic coefficients `[CFx, CFy, CFz, CMx, CMy, CMz]` of `body_aero` at angle of attack
`alpha` [rad], sideslip `beta` [rad], `wind_speed` [m/s] and rotation rate `body_aero.omega`
[rad/s], and their derivatives: `dalpha` and `dbeta` [1/rad], and `dp`, `dq`, `dr` with
respect to the body rates about x, y and z as p̂ = pb/2V, q̂ = q c_ref/2V and r̂ = rb/2V,
with b the wing span, c_ref `body_aero.c_ref` and V `wind_speed`. Moments are about, and
the body turns about, `solver.reference_point`. `kwargs` go to [`linearize`](@ref), which
leaves `body_aero` at this inflow.

Returns `(coeffs, dalpha, dbeta, dp, dq, dr, converged)`.
"""
function stability_derivatives(solver::Solver, body_aero::BodyAerodynamics, alpha, beta,
        wind_speed; kwargs...)
    va = apparent_wind(alpha, beta, wind_speed)
    set_va!(body_aero, va, body_aero.omega; reference_point=solver.reference_point)
    jac, results, converged = linearize(solver, body_aero, [va; body_aero.omega];
        theta_idxs=nothing, va_idxs=1:3, omega_idxs=4:6, aero_coeffs=true, kwargs...)
    dva_dalpha = ForwardDiff.derivative(
        angle -> apparent_wind(angle, beta, wind_speed), alpha)
    dva_dbeta = ForwardDiff.derivative(
        angle -> apparent_wind(alpha, angle, wind_speed), beta)
    coeff_jac = jac[1:6, :]
    span = body_aero.wings[1].span
    return (coeffs=results[1:6], dalpha=coeff_jac[:, 1:3] * dva_dalpha,
        dbeta=coeff_jac[:, 1:3] * dva_dbeta,
        dp=coeff_jac[:, 4] * 2wind_speed / span,
        dq=coeff_jac[:, 5] * 2wind_speed / body_aero.c_ref,
        dr=coeff_jac[:, 6] * 2wind_speed / span, converged)
end

"""
    trim_angle(solver, body_aero, beta, wind_speed; alpha_range=deg2rad.(-5:2:15),
               alpha_tol=1e-5, kwargs...)

Angles of attack [rad] at which `CMy` of `body_aero` about `solver.reference_point` changes
sign between neighbouring entries of `alpha_range`, bisected to `alpha_tol` [rad], at
sideslip `beta` [rad] and `wind_speed` [m/s]. Returns one `(alpha, dCMy_dalpha)` per trim,
the slope [1/rad] from [`stability_derivatives`](@ref) with `kwargs`; a trim is statically
stable where `dCMy_dalpha < 0`.
"""
function trim_angle(solver::Solver, body_aero::BodyAerodynamics, beta, wind_speed;
        alpha_range=deg2rad.(-5:2:15), alpha_tol=1e-5, kwargs...)
    is_nose_down =
        alpha -> pitch_moment_coeff(solver, body_aero, alpha, beta, wind_speed) < 0
    nose_down = is_nose_down.(alpha_range)
    trims = @NamedTuple{alpha::Float64, dCMy_dalpha::Float64}[]
    for i in 1:length(alpha_range)-1
        nose_down[i] == nose_down[i+1] && continue
        alpha = bisect_sign_change(is_nose_down, alpha_range[i], alpha_range[i+1],
            alpha_tol)
        derivatives = stability_derivatives(solver, body_aero, alpha, beta, wind_speed;
            kwargs...)
        push!(trims, (alpha=alpha, dCMy_dalpha=derivatives.dalpha[5]))
    end
    return trims
end

"""
    pitch_moment_coeff(solver, body_aero, alpha, beta, wind_speed)

`CMy` of `body_aero` solved at angle of attack `alpha` [rad], sideslip `beta` [rad] and
`wind_speed` [m/s], at the rotation rate `body_aero.omega` about `solver.reference_point`.
"""
function pitch_moment_coeff(solver, body_aero, alpha, beta, wind_speed)
    set_va!(body_aero, apparent_wind(alpha, beta, wind_speed), body_aero.omega;
        reference_point=solver.reference_point)
    return solve!(solver, body_aero).moment_coeffs[2]
end

"""
    bisect_sign_change(predicate, low, high, tol)

Bisect `[low, high]`, across which the boolean `predicate` flips, to a width of `tol` and
return the midpoint.
"""
function bisect_sign_change(predicate, low, high, tol)
    predicate_low = predicate(low)
    while high - low > tol
        middle = (low + high) / 2
        if predicate(middle) == predicate_low
            low = middle
        else
            high = middle
        end
    end
    return (low + high) / 2
end
