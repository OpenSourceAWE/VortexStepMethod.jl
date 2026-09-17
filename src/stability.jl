"""
    stability_derivatives(solver, body_aero, alpha, beta, wind_speed; kwargs...)

Aerodynamic coefficients `[CFx, CFy, CFz, CMx, CMy, CMz]` of `body_aero` at angle of attack
`alpha` [rad], sideslip `beta` [rad] and `wind_speed` [m/s], and their derivatives with
respect to `alpha` and `beta` [1/rad], at the rotation rate `body_aero.omega` and with
moments about `solver.reference_point`. `kwargs` go to [`linearize`](@ref), which leaves
`body_aero` at this inflow.

Returns `(coeffs, dalpha, dbeta, converged)`.
"""
function stability_derivatives(solver::Solver, body_aero::BodyAerodynamics, alpha, beta,
        wind_speed; kwargs...)
    va_vec = apparent_wind(alpha, beta, wind_speed)
    jac, results, converged = linearize(solver, body_aero, va_vec;
        theta_idxs=nothing, va_idxs=1:3, aero_coeffs=true, kwargs...)
    dva_dalpha = ForwardDiff.derivative(
        angle -> apparent_wind(angle, beta, wind_speed), alpha)
    dva_dbeta = ForwardDiff.derivative(
        angle -> apparent_wind(alpha, angle, wind_speed), beta)
    coeff_jac = jac[1:6, :]
    return (coeffs=results[1:6], dalpha=coeff_jac * dva_dalpha,
        dbeta=coeff_jac * dva_dbeta, converged)
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
`wind_speed` [m/s], at the rotation rate `body_aero.omega`.
"""
function pitch_moment_coeff(solver, body_aero, alpha, beta, wind_speed)
    set_va!(body_aero, apparent_wind(alpha, beta, wind_speed), body_aero.omega)
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
