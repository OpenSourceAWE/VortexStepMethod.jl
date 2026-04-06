using Test
using LinearAlgebra
using VortexStepMethod
using VortexStepMethod: Wing, Section, add_section!, refine!

"""
    build_flat_wing(n_panels, n_ribs; billowing_percentage=0.0)

Build a simple flat rectangular wing with `n_ribs` equally spaced rib
sections along the y-axis.  Chord runs from LE at x=0 to TE at x=-1.
"""
function build_flat_wing(n_panels, n_ribs;
                         billowing_percentage=0.0)
    dist = billowing_percentage > 0 ? BILLOWING : SPLIT_PROVIDED
    wing = Wing(n_panels;
        spanwise_distribution=dist,
        billowing_percentage=billowing_percentage)
    span = 10.0
    for i in 1:n_ribs
        y = span * (i - 1) / (n_ribs - 1)
        add_section!(wing,
            [0.0, y, 0.0],   # LE
            [-1.0, y, 0.0],  # TE
            INVISCID)
    end
    refine!(wing)
    return wing
end

"""
    te_arc_length_between_ribs(wing, rib_left, rib_right)

Sum up the TE segment lengths of refined sections between two
unrefined rib indices (inclusive of endpoints).
"""
function te_arc_length_between_ribs(wing, rib_left, rib_right)
    # Collect refined section indices belonging to this rib pair.
    # Rib sections sit at the boundaries; intermediate sections
    # are the ones with mapping equal to rib_left.
    sections = wing.refined_sections
    n = length(sections)
    # Find refined indices that map to this rib pair
    # The first refined section for rib_left is the rib itself,
    # then intermediate sections, then the next rib.
    first_idx = nothing
    last_idx = nothing
    for i in 1:n
        sec = sections[i]
        # Check if this section's LE matches the left rib
        if isnothing(first_idx)
            if isapprox(sec.LE_point,
                    wing.unrefined_sections[rib_left].LE_point;
                    atol=1e-10)
                first_idx = i
            end
        end
        if !isnothing(first_idx) && isnothing(last_idx)
            if isapprox(sec.LE_point,
                    wing.unrefined_sections[rib_right].LE_point;
                    atol=1e-10)
                last_idx = i
                break
            end
        end
    end
    @assert !isnothing(first_idx) && !isnothing(last_idx)
    arc = 0.0
    for i in first_idx:(last_idx - 1)
        arc += norm(sections[i + 1].TE_point -
                    sections[i].TE_point)
    end
    straight = norm(sections[last_idx].TE_point -
                    sections[first_idx].TE_point)
    return (; arc, straight)
end

@testset "Billowing" begin
    n_ribs = 5   # 4 rib pairs
    n_panels = 4 * 8  # 8 panels per rib pair

    @testset "Zero percentage = flat" begin
        wing = build_flat_wing(n_panels, n_ribs;
            billowing_percentage=0.0)
        for pair in 1:(n_ribs - 1)
            res = te_arc_length_between_ribs(wing, pair, pair + 1)
            @test res.arc ≈ res.straight atol=1e-10
        end
    end

    @testset "Rib endpoints preserved" begin
        wing = build_flat_wing(n_panels, n_ribs;
            billowing_percentage=10.0)
        for i in 1:n_ribs
            y_target = 10.0 * (i - 1) / (n_ribs - 1)
            # Find the refined section matching this rib
            found = false
            for sec in wing.refined_sections
                if isapprox(sec.LE_point,
                        wing.unrefined_sections[i].LE_point;
                        atol=1e-10)
                    @test isapprox(sec.TE_point,
                        wing.unrefined_sections[i].TE_point;
                        atol=1e-6)
                    found = true
                    break
                end
            end
            @test found
        end
    end

    @testset "Billowing sags outward" begin
        wing = build_flat_wing(n_panels, n_ribs;
            billowing_percentage=10.0)
        # For this flat test wing the negative rotation around
        # y_hat=[0,-1,0] produces -z sag in global coords.
        # On a real kite this maps to +z in the local airfoil frame.
        for pair in 1:(n_ribs - 1)
            le_left = wing.unrefined_sections[pair].LE_point
            le_right = wing.unrefined_sections[pair + 1].LE_point
            mid_y = (le_left[2] + le_right[2]) / 2
            best = nothing
            best_dist = Inf
            for sec in wing.refined_sections
                dist = abs(sec.LE_point[2] - mid_y)
                if dist < best_dist
                    best_dist = dist
                    best = sec
                end
            end
            @test best.TE_point[3] < 0
        end
    end

    @testset "Chord length preserved ($pct%)" for pct in
            [5.0, 10.0, 15.0]
        wing = build_flat_wing(n_panels, n_ribs;
            billowing_percentage=pct)
        for sec in wing.refined_sections
            chord = norm(sec.TE_point - sec.LE_point)
            # All ribs have chord=1.0, so interpolated chord=1.0
            @test isapprox(chord, 1.0; atol=1e-10)
        end
    end

    @testset "TE arc length matches percentage ($pct%)" for pct in
            [1.0, 5.0, 10.0, 15.0]
        wing = build_flat_wing(n_panels, n_ribs;
            billowing_percentage=pct)
        for pair in 1:(n_ribs - 1)
            res = te_arc_length_between_ribs(
                wing, pair, pair + 1)
            measured = (res.arc - res.straight) / res.arc * 100
            diff = measured - pct
            @info "Pair $pair: measured=$(round(measured; digits=3))%" *
                  " target=$(pct)% diff=$(round(diff; digits=3))%"
            @test isapprox(measured, pct; rtol=0.01)
        end
    end

    @testset "apply_billowing_to_pair! edge cases" begin
        using VortexStepMethod: apply_billowing_to_pair!, MVec3
        wing = build_flat_wing(n_panels, n_ribs;
            billowing_percentage=0.0)
        y_hat = MVec3(0.0, -1.0, 0.0)
        le_ref = MVec3(0.0, 2.5, 0.0)
        te_l = MVec3(-1.0, 0.0, 0.0)
        te_r = MVec3(-1.0, 2.5, 0.0)
        orig = [MVec3(s.TE_point)
                for s in wing.refined_sections]

        # Zero percentage is a no-op
        apply_billowing_to_pair!(
            wing.refined_sections, 2, 8,
            y_hat, 2.5, le_ref, te_l, te_r, 0.0)
        for (i, sec) in enumerate(wing.refined_sections)
            @test sec.TE_point == orig[i]
        end

        # Coincident ribs (straight ≈ 0) is a no-op
        apply_billowing_to_pair!(
            wing.refined_sections, 2, 8,
            y_hat, 2.5, le_ref, te_l, te_l, 10.0)
        for (i, sec) in enumerate(wing.refined_sections)
            @test sec.TE_point == orig[i]
        end
    end
end
