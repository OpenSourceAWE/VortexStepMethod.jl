"""
    KulfanBasis(; n_stations=60, n_weights=8)

The fixed CST basis a shape deformation is projected onto: chord stations `x` in
`[0, 1]` and the pseudoinverse `projection` of the class-times-Bernstein matrix
`C(x)·B(x)` those stations span.

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
    "`n_weights × length(x)` pseudoinverse of `C(x)·B(x)`, mapping deflection to weights."
    projection::Matrix{Float64}
    "Number of CST weights per surface, matching the airfoil being deformed."
    n_weights::Int
end

function KulfanBasis(; n_stations::Int=60, n_weights::Int=8)
    n_stations > n_weights || throw(ArgumentError(
        "KulfanBasis needs more stations than weights; got $n_stations and $n_weights."))
    theta = range(0, pi, n_stations)
    x = @. (1 - cos(theta)) / 2
    shape = class_function(x) .* bernstein_basis(x, n_weights - 1)
    return KulfanBasis(collect(x), pinv(shape), n_weights)
end

"""
    deform_kulfan(basis, base, upper_deflection, lower_deflection) -> KulfanParameters
    deform_kulfan(basis, base, camber) -> KulfanParameters

Add a surface deflection to a fixed set of Kulfan parameters, both deflections sampled
on `basis.x` and normalized by chord. The three-argument form applies one camber
deflection to both surfaces, leaving the thickness distribution untouched; the
four-argument form deforms the surfaces independently, which is what a double-skin
membrane with its own upper and lower control points needs.

The leading-edge weight and the trailing-edge thickness are carried over unchanged:
they are the two shape freedoms a chord-referenced deflection cannot resolve.
"""
function deform_kulfan(basis::KulfanBasis, base::KulfanParameters,
                       upper_deflection::AbstractVector,
                       lower_deflection::AbstractVector)
    length(base.upper_weights) == basis.n_weights || throw(ArgumentError(
        "KulfanBasis has $(basis.n_weights) weights, airfoil has " *
        "$(length(base.upper_weights))."))
    (length(upper_deflection) == length(basis.x) &&
     length(lower_deflection) == length(basis.x)) || throw(ArgumentError(
        "Deflections must be sampled on the basis' $(length(basis.x)) stations."))
    return KulfanParameters(base.upper_weights .+ basis.projection * upper_deflection,
                            base.lower_weights .+ basis.projection * lower_deflection,
                            base.leading_edge_weight, base.TE_thickness)
end

deform_kulfan(basis::KulfanBasis, base::KulfanParameters,
              camber::AbstractVector) = deform_kulfan(basis, base, camber, camber)

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
