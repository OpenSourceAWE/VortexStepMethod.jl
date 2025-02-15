import numpy as np
import pytest
from VSM.WingGeometry import Wing, Section  # Replace with your actual module names


def test_wing_initialization():
    example_wing = Wing(n_panels=10)
    assert example_wing.n_panels == 10
    assert example_wing.spanwise_panel_distribution == "linear"
    np.testing.assert_array_equal(example_wing.spanwise_direction, np.array([0, 1, 0]))
    assert isinstance(example_wing.sections, list)
    assert len(example_wing.sections) == 0  # Initially, sections list should be empty


def test_add_section():
    example_wing = Wing(n_panels=10)
    # Test adding a section
    example_wing.add_section(np.array([0, 0, 0]), np.array([-1, 0, 0]), "inviscid")
    assert len(example_wing.sections) == 1

    # Test if the section was added correctly
    section = example_wing.sections[0]
    np.testing.assert_array_equal(section.LE_point, np.array([0, 0, 0]))
    np.testing.assert_array_equal(section.TE_point, np.array([-1, 0, 0]))
    assert section.aero_input == "inviscid"


def test_robustness_left_to_right():
    example_wing = Wing(n_panels=10)
    # Test correct order
    example_wing.add_section(np.array([0, 1, 0]), np.array([0, 1, 0]), ["inviscid"])
    example_wing.add_section(np.array([0, -1, 0]), np.array([0, -1, 0]), ["inviscid"])
    example_wing.add_section(np.array([0, -1.5, 0]), np.array([0, -1.5, 0]), ["inviscid"])
    example_wing.refine_aerodynamic_mesh()

    example_wing_1 = Wing(n_panels=10)
    # Test right to left order
    example_wing_1.add_section(np.array([0, -1.5, 0]), np.array([0, -1.5, 0]), ["inviscid"])
    example_wing_1.add_section(np.array([0, -1, 0]), np.array([0, -1, 0]), ["inviscid"])
    example_wing_1.add_section(np.array([0, 1, 0]), np.array([0, 1, 0]), ["inviscid"])
    example_wing_1.refine_aerodynamic_mesh()

    example_wing_2 = Wing(n_panels=10)
    # Test random order
    example_wing_2.add_section(np.array([0, 1, 0]), np.array([0, 1, 0]), ["inviscid"])
    example_wing_2.add_section(np.array([0, -1.5, 0]), np.array([0, -1.5, 0]), ["inviscid"])
    example_wing_2.add_section(np.array([0, -1, 0]), np.array([0, -1, 0]), ["inviscid"])
    example_wing_2.refine_aerodynamic_mesh()

    for i in range(len(example_wing.sections)):
        np.testing.assert_array_equal(
            example_wing.sections[i].LE_point, example_wing_1.sections[i].LE_point
        )
        np.testing.assert_array_equal(
            example_wing.sections[i].TE_point, example_wing_1.sections[i].TE_point
        )
        np.testing.assert_array_equal(
            example_wing.sections[i].LE_point, example_wing_2.sections[i].LE_point
        )
        np.testing.assert_array_equal(
            example_wing.sections[i].TE_point, example_wing_2.sections[i].TE_point
        )


def test_refine_aerodynamic_mesh():
    n_panels = 4
    span = 20

    ## Test linear distribution
    wing = Wing(n_panels, spanwise_panel_distribution="linear")
    wing.add_section([0, span / 2, 0], [-1, span / 2, 0], ["inviscid"])
    wing.add_section([0, -span / 2, 0], [-1, -span / 2, 0], ["inviscid"])
    sections = wing.refine_aerodynamic_mesh()

    # Each panel has 2 corresponding sections, so we expect n_panels + 1 sections
    assert len(sections) == wing.n_panels + 1

    # Numerical tests for analytical calculations
    # Assuming linear interpolation, check specific points
    for i in range(len(sections)):
        # Calculate expected points for linear interpolation
        expected_LE = np.array([0, span / 2 - i * span / n_panels, 0])
        expected_TE = np.array([-1, span / 2 - i * span / n_panels, 0])

        # Check interpolated points with adjusted tolerance
        np.testing.assert_allclose(sections[i].LE_point, expected_LE, rtol=1e-5)
        np.testing.assert_allclose(sections[i].TE_point, expected_TE, rtol=1e-4)

    ## Test cosine distribution
    wing = Wing(n_panels, spanwise_panel_distribution="cosine")
    wing.add_section([0, span / 2, 0], [-1, span / 2, 0], ["inviscid"])
    wing.add_section([0, -span / 2, 0], [-1, -span / 2, 0], ["inviscid"])
    sections = wing.refine_aerodynamic_mesh()
    print(sections[0].LE_point)
    # Each panel has 2 corresponding sections, so we expect n_panels + 1 sections
    assert len(sections) == wing.n_panels + 1

    # Numerical tests for cosine interpolation
    for i in range(len(sections)):
        # Cosine distribution formula for linear span
        theta = np.linspace(0, np.pi, n_panels + 1)
        expected_LE_y = span / 2 * np.cos(theta)
        expected_LE = np.array([0, expected_LE_y[i], 0])
        expected_TE = np.array([-1, expected_LE_y[i], 0])

        # Check interpolated points with adjusted tolerance
        np.testing.assert_allclose(sections[i].LE_point, expected_LE, atol=1e-8)
        np.testing.assert_allclose(sections[i].TE_point, expected_TE, atol=1e-8)


def test_refine_aerodynamic_mesh_1_panel():
    n_panels = 1
    span = 20

    wing = Wing(n_panels, spanwise_panel_distribution="linear")
    wing.add_section([0, span / 2, 0], [-1, span / 2, 0], ["inviscid"])
    wing.add_section([0, -span / 2, 0], [-1, -span / 2, 0], ["inviscid"])

    sections = wing.refine_aerodynamic_mesh()

    # Each panel has 2 corresponding sections, so we expect n_panels + 1 sections
    assert len(sections) == wing.n_panels + 1

    np.testing.assert_array_equal(sections[0].LE_point, np.array([0, span / 2, 0]))
    np.testing.assert_array_equal(sections[0].TE_point, np.array([-1, span / 2, 0]))


def test_refine_aeordynamic_mesh_2_panel():
    n_panels = 2
    span = 20

    wing = Wing(n_panels, spanwise_panel_distribution="linear")
    wing.add_section([0, span / 2, 0], [-1, span / 2, 0], ["inviscid"])
    wing.add_section([0, -span / 2, 0], [-1, -span / 2, 0], ["inviscid"])

    sections = wing.refine_aerodynamic_mesh()

    # Each panel has 2 corresponding sections, so we expect n_panels + 1 sections
    assert len(sections) == wing.n_panels + 1

    np.testing.assert_array_equal(sections[0].LE_point, np.array([0, span / 2, 0]))
    np.testing.assert_array_equal(sections[1].LE_point, np.array([0, 0, 0]))
    np.testing.assert_array_equal(sections[2].LE_point, np.array([0, -span / 2, 0]))


def test_refine_aeordynamic_mesh_more_sections_than_panels():
    n_panels = 2
    span = 20

    wing = Wing(n_panels, spanwise_panel_distribution="linear")
    wing.add_section([0, span / 2, 0], [-1, span / 2, 0], ["inviscid"])
    wing.add_section([0, span / 4, 0], [-1, span / 4, 0], ["inviscid"])
    wing.add_section([0, 0, 0], [-1, 0, 0], ["inviscid"])
    wing.add_section([0, -span / 4, 0], [-1, -span / 4, 0], ["inviscid"])
    wing.add_section([0, -span / 3, 0], [-1, -span / 3, 0], ["inviscid"])
    wing.add_section([0, -span / 2, 0], [-1, -span / 2, 0], ["inviscid"])

    sections = wing.refine_aerodynamic_mesh()

    # Each panel has 2 corresponding sections, so we expect n_panels + 1 sections
    assert len(sections) == wing.n_panels + 1

    for i in range(len(sections)):
        # Calculate expected points for linear interpolation
        expected_LE = np.array([0, span / 2 - i * span / n_panels, 0])
        expected_TE = np.array([-1, span / 2 - i * span / n_panels, 0])

        # Check interpolated points with adjusted tolerance
        np.testing.assert_allclose(sections[i].LE_point, expected_LE, rtol=1e-5)
        np.testing.assert_allclose(sections[i].TE_point, expected_TE, rtol=1e-4)


def test_refine_aerodynamic_mesh_for_symmetrical_wing():
    n_panels = 2
    span = 10  # Total span from -5 to 5

    wing = Wing(n_panels, spanwise_panel_distribution="linear")
    wing.add_section([0, 5, 0], [-1, 5, 0], ["inviscid"])
    wing.add_section([0, -5, 0], [-1, -5, 0], ["inviscid"])


    sections = wing.refine_aerodynamic_mesh()

    # Calculate expected quarter-chord points
    qc_start = np.array([0.25 * -1, 5, 0])
    qc_end = np.array([0.25 * -1, -5, 0])
    expected_qc_y = np.linspace(qc_start[1], qc_end[1], n_panels + 1)

    for i, section in enumerate(sections):
        # Calculate expected quarter-chord point
        expected_qc = np.array([0.25 * -1, expected_qc_y[i], 0])

        # Calculate expected chord vector
        chord_start = np.array([-1, 5, 0]) - np.array([0, 5, 0])
        chord_end = np.array([-1, -5, 0]) - np.array([0, -5, 0])
        t = (expected_qc_y[i] - qc_start[1]) / (qc_end[1] - qc_start[1])

        # Normalize chord vectors
        chord_start_norm = chord_start / np.linalg.norm(chord_start)
        chord_end_norm = chord_end / np.linalg.norm(chord_end)

        # Interpolate direction
        avg_direction = (1 - t) * chord_start_norm + t * chord_end_norm
        avg_direction = avg_direction / np.linalg.norm(avg_direction)

        # Interpolate length
        chord_start_length = np.linalg.norm(chord_start)
        chord_end_length = np.linalg.norm(chord_end)
        avg_length = (1 - t) * chord_start_length + t * chord_end_length

        expected_chord = avg_direction * avg_length

        # Calculate expected LE and TE points
        expected_LE = expected_qc - 0.25 * expected_chord
        expected_TE = expected_qc + 0.75 * expected_chord

        print(f"Section {i}:")
        print(f"  Expected LE: {expected_LE}")
        print(f"  Actual LE:   {section.LE_point}")
        print(f"  Expected TE: {expected_TE}")
        print(f"  Actual TE:   {section.TE_point}")
        print(f"  Expected chord: {expected_chord}")
        print(f"  Actual chord:   {section.TE_point - section.LE_point}")

        # Check interpolated points with adjusted tolerance
        np.testing.assert_allclose(section.LE_point, expected_LE, rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(section.TE_point, expected_TE, rtol=1e-5, atol=1e-5)

    # Additional checks
    assert len(sections) == n_panels + 1
    assert sections[0].LE_point[1] == pytest.approx(5, abs=1e-5)
    assert sections[-1].LE_point[1] == pytest.approx(-5, abs=1e-5)
    assert sections[0].TE_point[1] == pytest.approx(5, abs=1e-5)
    assert sections[-1].TE_point[1] == pytest.approx(-5, abs=1e-5)


def test_refine_aeordynamic_mesh_lei_airfoil_interpolation():
    n_panels = 4
    span = 20

    wing = Wing(n_panels, spanwise_panel_distribution="linear")
    wing.add_section(
        [0, span / 2, 0], [-1, span / 2, 0], ["lei_airfoil_breukels", [0, 0]]
    )
    wing.add_section([0, 0, 0], [-1, 0, 0], ["lei_airfoil_breukels", [2, 0.5]])
    wing.add_section(
        [0, -span / 2, 0], [-1, -span / 2, 0], ["lei_airfoil_breukels", [4, 1]]
    )

    sections = wing.refine_aerodynamic_mesh()

    # Each panel has 2 corresponding sections, so we expect n_panels + 1 sections
    assert len(sections) == wing.n_panels + 1

    expected_tube_diameter = np.linspace(0, 4, n_panels + 1)
    expected_chamber_height = np.linspace(0, 1, n_panels + 1)
    print(f"expected_tube_diameter: {expected_tube_diameter}")
    print(f"expected_chamber_height: {expected_chamber_height}")

    for i, section in enumerate(sections):
        # Calculate expected points for linear interpolation
        expected_LE = np.array([0, span / 2 - i * span / n_panels, 0])
        expected_TE = np.array([-1, span / 2 - i * span / n_panels, 0])

        # Check interpolated points with adjusted tolerance
        np.testing.assert_allclose(section.LE_point, expected_LE, rtol=1e-5)
        np.testing.assert_allclose(section.TE_point, expected_TE, rtol=1e-4)

        # Check if the airfoil data is correctly interpolated
        aero_input_i = section.aero_input
        assert np.array_equal(str(aero_input_i[0]), "lei_airfoil_breukels")

        tube_diam = aero_input_i[1][0]
        chamber_height = aero_input_i[1][1]
        assert np.isclose(tube_diam, expected_tube_diameter[i])
        assert np.isclose(chamber_height, expected_chamber_height[i])


def test_refine_mesh_by_splitting_provided_sections():
    # Create mock sections
    section1 = Section(
        LE_point=np.array([0, 0, 0]),
        TE_point=np.array([1, 0, 0]),
        aero_input=["inviscid"],
    )
    section2 = Section(
        LE_point=np.array([0, 1, 0]),
        TE_point=np.array([1, 1, 0]),
        aero_input=["inviscid"],
    )
    section3 = Section(
        LE_point=np.array([0, 2, 0]),
        TE_point=np.array([1, 2, 0]),
        aero_input=["inviscid"],
    )

    # Create mock Wing object
    wing = Wing(
        sections=[section1, section2, section3],
        n_panels=6,
        spanwise_panel_distribution="split_provided",
    )

    # Call the function
    new_sections = wing.refine_mesh_by_splitting_provided_sections()

    # Assert the correct number of sections
    assert (
        len(new_sections) - 1 == 6
    ), f"Expected 6 sections, but got {len(new_sections)}"

    # Check if the original sections are preserved
    assert new_sections[0] == section1, "First section should be preserved"
    assert new_sections[3] == section2, "Middle section should be preserved"
    assert new_sections[-1] == section3, "Last section should be preserved"

    # Check if new sections are created between original sections
    assert (
        new_sections[1].LE_point[1] > 0 and new_sections[1].LE_point[1] < 1
    ), "New section 1 should be between section 1 and 2"
    assert (
        new_sections[2].LE_point[1] > 0 and new_sections[2].LE_point[1] < 1
    ), "New section 2 should be between section 1 and 2"
    assert (
        new_sections[4].LE_point[1] > 1 and new_sections[4].LE_point[1] < 2
    ), "New section 4 should be between section 2 and 3"

    # Check if new sections have correct aero_input
    for section in new_sections:
        assert section.aero_input == [
            "inviscid"
        ], f"Expected aero_input ['inviscid'], but got {section.aero_input}"


if __name__ == "__main__":
    pytest.main()
