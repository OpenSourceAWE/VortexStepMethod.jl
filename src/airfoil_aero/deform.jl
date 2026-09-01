"""
    KulfanBasis(; n_stations=60, n_weights=8, ridge=1e-3)

The fixed CST basis a shape deformation is projected onto: chord stations `x` in
`[0, 1]` and the ridge-regularised inverse `projection` of the class-times-Bernstein
matrix `C(x)·B(x)` those stations span.

`ridge` is the Tikhonov weight, relative to the basis' own largest singular value, that
keeps the projection from answering a deflection it cannot represent with weights far
larger than the airfoil they correct. A plain pseudoinverse has no such bound: a kinked
deflection — a strut buckling is one — lands on the basis' weakest directions and comes
back amplified a hundredfold and alternating in sign, which is not an airfoil. The ridge
trades a small, measured under-response on deflections the basis *can* hold for a bounded
answer on the ones it cannot.

CST is linear in its weights, so a surface displacement is a matvec against this
constant matrix — never a refit. Refitting inside a loop is not an option: the Kulfan
fit is non-unique, so the same shape fitted twice returns weight vectors differing by
more than the deformation signal, which reaches the polars as frame-to-frame jitter.

`C(1) = 0`, so the parameterisation cannot put the trailing edge off the chord line.
Deflections must therefore be measured **against the deformed chord**, with the chord
rotation and stretch taken from the leading- and trailing-edge points; what is left is
the representable residual.
"""
struct KulfanBasis
    "Chord stations the deflection is sampled on, ascending in `[0, 1]`."
    x::Vector{Float64}
    "`n_weights × length(x)` regularised inverse of `C(x)·B(x)`, mapping deflection to weights."
    projection::Matrix{Float64}
    "Number of CST weights per surface, matching the airfoil being deformed."
    n_weights::Int
    "Tikhonov weight the projection was built with, relative to the basis' largest singular value."
    ridge::Float64
end

function KulfanBasis(; n_stations::Int=60, n_weights::Int=8, ridge::Real=1e-3)
    n_stations > n_weights || throw(ArgumentError(
        "KulfanBasis needs more stations than weights; got $n_stations and $n_weights."))
    ridge >= 0 || throw(ArgumentError("KulfanBasis ridge must not be negative."))
    theta = range(0, pi, n_stations)
    x = @. (1 - cos(theta)) / 2
    shape = class_function(x) .* bernstein_basis(x, n_weights - 1)
    scale = maximum(svdvals(shape))^2
    projection = iszero(ridge) ? pinv(shape) :
        (shape' * shape + ridge * scale * I) \ shape'
    return KulfanBasis(collect(x), projection, n_weights, Float64(ridge))
end

"""
    deform_kulfan(basis, base, upper_deflection, lower_deflection) -> KulfanParameters
    deform_kulfan(basis, base, camber) -> KulfanParameters

Add a surface deflection to a fixed set of Kulfan parameters, both deflections sampled
on `basis.x` and normalized by chord. The three-argument form applies one camber
deflection to both surfaces, leaving the thickness distribution untouched; the
four-argument form deforms the surfaces independently, which is what a double-skin
membrane with its own upper and lower control points needs.

Each deflection is first reduced to the part the basis can carry, see
[`chord_residual`](@ref): the straight line through its own endpoints is a chord
rotation and translation, which the basis cannot express and `pinv` answers with
runaway weights. Take that line from [`chord_line`](@ref) if the frame the deflection
was measured in has not already absorbed it.

The leading-edge weight and the trailing-edge thickness are carried over unchanged:
they are the two shape freedoms a chord-referenced deflection cannot resolve.
"""
function deform_kulfan(basis::KulfanBasis, base::KulfanParameters,
                       upper_deflection::AbstractVector,
                       lower_deflection::AbstractVector)
    out = KulfanParameters(copy(base.upper_weights), copy(base.lower_weights),
                           base.leading_edge_weight, base.TE_thickness)
    return deform_kulfan!(out, basis, base, upper_deflection, lower_deflection,
                          similar(basis.x))
end

"""
    deform_kulfan!(out, basis, base, upper_deflection, lower_deflection, residual)
        -> out

[`deform_kulfan`](@ref) written into the shape `out` already is, with `residual` as the
scratch [`chord_residual!`](@ref) fills — the form a live polar refreshes a panel's shape
with, since it neither allocates nor replaces the object the panel points at.
"""
function deform_kulfan!(out::KulfanParameters, basis::KulfanBasis,
                        base::KulfanParameters, upper_deflection::AbstractVector,
                        lower_deflection::AbstractVector, residual::AbstractVector)
    length(base.upper_weights) == basis.n_weights || throw(ArgumentError(
        "KulfanBasis has $(basis.n_weights) weights, airfoil has " *
        "$(length(base.upper_weights))."))
    (length(upper_deflection) == length(basis.x) &&
     length(lower_deflection) == length(basis.x)) || throw(ArgumentError(
        "Deflections must be sampled on the basis' $(length(basis.x)) stations."))
    out.leading_edge_weight = base.leading_edge_weight
    out.TE_thickness = base.TE_thickness
    for (weights, base_weights, deflection) in
            ((out.upper_weights, base.upper_weights, upper_deflection),
             (out.lower_weights, base.lower_weights, lower_deflection))
        chord_residual!(residual, basis, deflection)
        mul!(weights, basis.projection, residual)
        weights .+= base_weights
    end
    return out
end

"""
    chord_residual(basis, deflection) -> Vector{Float64}

The part of a deflection the CST basis can represent. `C(x)` vanishes at both chord
ends, so the basis can put neither the leading nor the trailing edge off the chord
line: a deflection that does not end at zero is asking for a chord **rotation and
translation**, not a camber change. Projecting it anyway is not merely inexact, it is
unstable — `pinv` answers an unrepresentable end displacement with weights an order of
magnitude past the ones it is correcting, and the airfoil that comes back is not one.

Removes the straight line through the deflection's own endpoints, which
[`chord_line`](@ref) returns, and gives back what is left.
"""
chord_residual(basis::KulfanBasis, deflection::AbstractVector) =
    chord_residual!(similar(basis.x), basis, deflection)

"""
    chord_residual!(residual, basis, deflection) -> residual

[`chord_residual`](@ref) written into storage the caller owns.
"""
function chord_residual!(residual::AbstractVector, basis::KulfanBasis,
                         deflection::AbstractVector)
    offset, slope = chord_line(basis, deflection)
    residual .= deflection .- (offset .+ slope .* basis.x)
    return residual
end

"""
    chord_line(basis, deflection) -> (offset, slope)

The straight line through a deflection's own endpoints, both over chord: `offset` moves
the leading edge and `slope` is the chord rotation `atan(slope)` [rad] the caller owes
its angle of attack, when the frame it measured the deflection in has not already
absorbed it. It is the part [`chord_residual`](@ref) removes.
"""
function chord_line(basis::KulfanBasis, deflection::AbstractVector)
    length(deflection) == length(basis.x) || throw(ArgumentError(
        "Deflection must be sampled on the basis' $(length(basis.x)) stations."))
    span = basis.x[end] - basis.x[1]
    abs(span) < eps() && return (deflection[1], 0.0)
    slope = (deflection[end] - deflection[1]) / span
    return (deflection[1] - slope * basis.x[1], slope)
end

deform_kulfan(basis::KulfanBasis, base::KulfanParameters,
              camber::AbstractVector) = deform_kulfan(basis, base, camber, camber)

deform_kulfan!(out::KulfanParameters, basis::KulfanBasis, base::KulfanParameters,
               camber::AbstractVector, residual::AbstractVector) =
    deform_kulfan!(out, basis, base, camber, camber, residual)

"""
    control_point_deflection(basis, fractions, deflections) -> Vector{Float64}

Resample a deflection given at arbitrary chord `fractions` onto `basis.x`, by linear
interpolation, held flat outside the sampled range. This is the generic bridge from
control points — beam nodes now, membrane nodes later — to the CST basis: the caller
only has to say where along the chord each point sits and how far off the chord line it
moved, both normalized by the local chord.

`fractions` need not be sorted, but must not repeat a station.
"""
function control_point_deflection(basis::KulfanBasis, fractions::AbstractVector,
                                  deflections::AbstractVector)
    length(fractions) == length(deflections) || throw(ArgumentError(
        "control_point_deflection: $(length(fractions)) fractions for " *
        "$(length(deflections)) deflections."))
    order = sortperm(collect(float.(fractions)))
    knots = collect(float.(fractions))[order]
    values = collect(float.(deflections))[order]
    if length(knots) < 2
        return fill(isempty(values) ? 0.0 : values[1], length(basis.x))
    end
    allunique(knots) || throw(ArgumentError(
        "control_point_deflection: repeated chord fraction in $knots."))
    interp = linear_interpolation(knots, values; extrapolation_bc=Flat())
    return [interp(xi) for xi in basis.x]
end
