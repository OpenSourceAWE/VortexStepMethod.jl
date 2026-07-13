"""
    CpData

Raw chord-slice surface-pressure samples for one 2D section (no interpolants) — the
CSV/section-level counterpart of [`CpPolar`](@ref). Sampled at `n_chord` **vertical
chord slices** `chord_x` (x/c), over angle of attack `alpha_range` and trailing-edge
deflection `delta_range` (both radians).

Fields:
- `n_chord`, `chord_x`: number and x/c of the chord slices.
- `alpha_range`, `delta_range`: grid axes (radians).
- `cp_upper`, `cp_lower`: samples, `n_chord × n_alpha × n_delta`.
"""
struct CpData
    n_chord::Int
    chord_x::Vector{Float64}
    alpha_range::Vector{Float64}
    delta_range::Vector{Float64}
    cp_upper::Array{Float64,3}
    cp_lower::Array{Float64,3}
end

"""
    CpPolar

Interpolatable surface-pressure table for a 2D section: a [`CpData`](@ref) plus one
`cp(alpha, delta)` interpolant per chord slice for each surface. `ΔCp = Cp_lower -
Cp_upper` is the chordwise load. Build one with `CpPolar(::CpData)`.

Fields:
- `data`: the underlying [`CpData`](@ref).
- `cp_upper_interp`, `cp_lower_interp`: one `cp(alpha, delta)` interpolant per slice.
"""
struct CpPolar{I}
    data::CpData
    cp_upper_interp::Vector{I}
    cp_lower_interp::Vector{I}
end
