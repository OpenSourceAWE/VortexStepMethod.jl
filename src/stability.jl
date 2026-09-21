"""
    stability_derivatives(solver, body_aero, alpha, beta, va; kwargs...)

Aerodynamic coefficients `[CFx, CFy, CFz, CMx, CMy, CMz]` of `body_aero` at angle of attack
`alpha` [rad], sideslip `beta` [rad], apparent wind speed `va` [m/s] and rotation rate
`body_aero.omega` [rad/s], and their derivatives: `dalpha` and `dbeta` [1/rad], and `dp`,
`dq`, `dr` with respect to p̂ = pb/2V, q̂ = q c_ref/2V and r̂ = rb/2V, where p, q, r are the
rates about body x, y, z, b the wing span, c_ref `body_aero.c_ref` and V is `va`, not the
area-weighted panel inflow speed the coefficients are normalised by. Moments are about,
and the body turns about, `solver.reference_point`. `kwargs` go to [`linearize`](@ref).
Leaves `body_aero` at this inflow with `body_aero.reference_point = solver.reference_point`.

Returns `(coeffs, dalpha, dbeta, dp, dq, dr, converged)`.
"""
function stability_derivatives(solver::Solver, body_aero::BodyAerodynamics, alpha, beta,
        va; kwargs...)
    va_vec = apparent_wind(alpha, beta, va)
    set_va!(body_aero, va_vec, body_aero.omega; reference_point=solver.reference_point)
    jac, results, converged = linearize(solver, body_aero, [va_vec; body_aero.omega];
        theta_idxs=nothing, va_vec_idxs=1:3, omega_idxs=4:6, aero_coeffs=true, kwargs...)
    dva_dalpha = ForwardDiff.derivative(
        angle -> apparent_wind(angle, beta, va), alpha)
    dva_dbeta = ForwardDiff.derivative(
        angle -> apparent_wind(alpha, angle, va), beta)
    coeff_jac = jac[1:6, :]
    span = body_aero.wings[1].span
    return (coeffs=results[1:6], dalpha=coeff_jac[:, 1:3] * dva_dalpha,
        dbeta=coeff_jac[:, 1:3] * dva_dbeta,
        dp=coeff_jac[:, 4] * 2va / span,
        dq=coeff_jac[:, 5] * 2va / body_aero.c_ref,
        dr=coeff_jac[:, 6] * 2va / span, converged)
end

"""
    trim_angle(solver, body_aero, beta, va; alpha_range=deg2rad.(-5:2:15),
               alpha_tol=1e-5, backend=AutoForwardDiff())

Angles of attack [rad] at which `CMy` of `body_aero` about `solver.reference_point` changes
sign between neighbouring entries of `alpha_range`, bisected to `alpha_tol` [rad], at
sideslip `beta` [rad] and apparent wind speed `va` [m/s]. Returns one
`(alpha, dCMy_dalpha)` per trim, the slope [1/rad] from [`stability_derivatives`](@ref)
with `backend`; a trim is statically stable where `dCMy_dalpha < 0`. Throws a
[`SolveFailure`](@ref) if a solve misses the solver's tolerances. Leaves
`body_aero.reference_point = solver.reference_point`.
"""
function trim_angle(solver::Solver, body_aero::BodyAerodynamics, beta, va;
        alpha_range=deg2rad.(-5:2:15), alpha_tol=1e-5, backend=AutoForwardDiff())
    is_nose_down = alpha -> nose_down(solver, body_aero, alpha, beta, va)
    nose_down_range = is_nose_down.(alpha_range)
    trims = @NamedTuple{alpha::Float64, dCMy_dalpha::Float64}[]
    for i in 1:length(alpha_range)-1
        nose_down_range[i] == nose_down_range[i+1] && continue
        alpha = bisect_sign_change(is_nose_down, alpha_range[i], alpha_range[i+1],
            alpha_tol)
        derivatives = stability_derivatives(solver, body_aero, alpha, beta, va;
            backend, throw_on_fail=true)
        push!(trims, (alpha=alpha, dCMy_dalpha=derivatives.dalpha[5]))
    end
    return trims
end

"""
    coeffs_at_angles(solver, body_aero, alpha, beta, va)

Aerodynamic coefficients `[CFx, CFy, CFz, CMx, CMy, CMz]` of `body_aero` solved at angle of
attack `alpha` [rad], sideslip `beta` [rad] and apparent wind speed `va` [m/s], at the
rotation rate `body_aero.omega`, which it turns about and stores as
`body_aero.reference_point`. Throws a [`SolveFailure`](@ref) if the solve misses the
solver's tolerances.
"""
function coeffs_at_angles(solver, body_aero, alpha, beta, va)
    set_va!(body_aero, apparent_wind(alpha, beta, va), body_aero.omega;
        reference_point=solver.reference_point)
    sol = solve!(solver, body_aero; throw_on_fail=true)
    return [sol.force_coeffs; sol.moment_coeffs]
end

"""
    nose_down(solver, body_aero, alpha, beta, va)

Whether `CMy` from [`coeffs_at_angles`](@ref) is negative.
"""
nose_down(solver, body_aero, alpha, beta, va) =
    coeffs_at_angles(solver, body_aero, alpha, beta, va)[5] < 0

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
