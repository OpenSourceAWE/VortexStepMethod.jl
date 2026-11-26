using Test, LinearAlgebra, VortexStepMethod

@testset "kite_sections builds the kite its KiteShape describes" begin
    shape = KiteShape(projected_area=12.0, aspect_ratio=4.0, arc_angle=70.0,
                      arc_exponent=1.5, taper=0.4, sweep=20.0, twist=5.0)
    flat = KiteShape(projected_area=12.0, aspect_ratio=4.0, arc_angle=70.0,
                     arc_exponent=1.5, taper=0.4, sweep=20.0, twist=0.0)
    sections = kite_sections(shape; n_sections=201)
    tip, centre = 1, 101
    span = sqrt(shape.aspect_ratio * shape.projected_area)
    chord(i) = norm(sections.te[i] - sections.le[i])

    @testset "projected span and area" begin
        @test sections.le[tip][2] - sections.le[end][2] ≈ span rtol = 1e-6
        wing = Wing(200; spanwise_distribution=UNCHANGED)
        flat_sections = kite_sections(flat; n_sections=201)
        foreach((le, te) -> add_section!(wing, le, te, INVISCID), flat_sections.le,
                flat_sections.te)
        refine!(wing)
        @test calculate_projected_area(wing) ≈ shape.projected_area rtol = 1e-3
    end

    @testset "arc, taper, sweep and twist at the tip" begin
        arc = sections.le[tip] - sections.le[tip + 1]
        midway = (sections.s[tip] + sections.s[tip + 1]) / 2
        @test rad2deg(atan(-arc[3], arc[2])) ≈
              shape.arc_angle * midway^shape.arc_exponent atol = 1e-3
        @test chord(tip) / chord(centre) ≈ shape.taper rtol = 1e-9
        @test sections.le[tip][1] - sections.le[centre][1] ≈
              tan(deg2rad(shape.sweep)) * span / 2 rtol = 1e-9
        normal = normalize(cross(MVec3(1, 0, 0), arc))
        @test rad2deg(asin(dot(sections.te[tip] - sections.le[tip], normal) /
                           chord(tip))) ≈ shape.twist atol = 0.05
        @test sections.te[centre] - sections.le[centre] ≈ [chord(centre), 0, 0]
    end

    @testset "mirror symmetric about the centre plane" begin
        mirror(point) = MVec3(point[1], -point[2], point[3])
        @test all(mirror.(sections.le) .≈ reverse(sections.le))
        @test all(mirror.(sections.te) .≈ reverse(sections.te))
    end
end
