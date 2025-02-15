import pytest
import numpy as np
import logging
from VSM.Panel import Panel  # Assuming the Panel class is in a file named Panel.py

import os
import sys

# Go back to root folder
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, root_path)
import tests.utils as test_utils
import tests.thesis_functions_oriol_cayon as thesis_functions


# Mock Section class for testing
class MockSection:
    def __init__(self, LE_point, TE_point, aero_input):
        self.LE_point = np.array(LE_point)
        self.TE_point = np.array(TE_point)
        self.aero_input = aero_input


def create_panel(section1, section2):
    section = {
        "p1": section1.LE_point,
        "p2": section2.LE_point,
        "p3": section2.TE_point,
        "p4": section1.TE_point,
    }
    bound_1 = section["p1"] * 3 / 4 + section["p4"] * 1 / 4
    bound_2 = section["p2"] * 3 / 4 + section["p3"] * 1 / 4

    mid_LE_point = section2.LE_point + 0.5 * (section1.LE_point - section2.LE_point)
    mid_TE_point = section2.TE_point + 0.5 * (section1.TE_point - section2.TE_point)
    mid_LE_vector = mid_TE_point - mid_LE_point
    aerodynamic_center = bound_1 + 0.5 * (bound_2 - bound_1)
    control_point = aerodynamic_center + 0.5 * mid_LE_vector

    LLpoint = aerodynamic_center
    VSMpoint = control_point
    x_airf = np.cross(VSMpoint - LLpoint, section["p2"] - section["p1"])
    x_airf = x_airf / np.linalg.norm(x_airf)

    # TANGENTIAL y_airf defined parallel to the chord-line, from LE-to-TE
    y_airf = VSMpoint - LLpoint
    y_airf = y_airf / np.linalg.norm(y_airf)

    # SPAN z_airf along the LE, in plane (towards left tip, along span) from the airfoil perspective
    z_airf = bound_2 - bound_1
    z_airf = z_airf / np.linalg.norm(z_airf)

    return Panel(
        section1,
        section2,
        aerodynamic_center,
        control_point,
        bound_1,
        bound_2,
        x_airf,
        y_airf,
        z_airf,
    )


@pytest.fixture
def sample_panel():
    section1 = MockSection([0, 0, 0], [1, 0, 0], ["inviscid"])
    section2 = MockSection([0, 10, 0], [1, 10, 0], ["inviscid"])
    return create_panel(section1, section2)


def test_panel_initialization(sample_panel):
    assert isinstance(sample_panel, Panel)


def test_panel_te_le_points(sample_panel):
    assert np.allclose(sample_panel.TE_point_1, [1, 0, 0])
    assert np.allclose(sample_panel.TE_point_2, [1, 10, 0])
    assert np.allclose(sample_panel.LE_point_1, [0, 0, 0])
    assert np.allclose(sample_panel.LE_point_2, [0, 10, 0])


def test_panel_corner_points(sample_panel):
    expected_corner_points = np.array(
        [
            [0, 0, 0],  # LE_point_1
            [1, 0, 0],  # TE_point_1
            [1, 10, 0],  # TE_point_2
            [0, 10, 0],  # LE_point_2
        ]
    )
    assert np.allclose(sample_panel.corner_points, expected_corner_points)


def test_panel_chord(sample_panel):
    rib_1 = sample_panel.corner_points[1] - sample_panel.corner_points[0]
    norm_rib_1 = np.linalg.norm(rib_1)
    rib_2 = sample_panel.corner_points[3] - sample_panel.corner_points[2]
    norm_rib_2 = np.linalg.norm(rib_2)
    chord_sample_panel = (norm_rib_1 + norm_rib_2) / 2
    assert np.isclose(sample_panel.chord, chord_sample_panel)


def test_va_initialization(sample_panel):
    assert sample_panel.va is None


# testing polar_data_input_option
def test_polar_data_input():
    # Generate mock polar data, using inviscid standards
    aoa = np.arange(-180, 180, 1)
    airfoil_data = np.empty((len(aoa), 4))
    for j, alpha in enumerate(aoa):
        cl, cd, cm = 2 * np.pi * np.deg2rad(alpha), 0.05, 0.01
        airfoil_data[j, 0] = np.deg2rad(alpha)
        airfoil_data[j, 1] = cl
        airfoil_data[j, 2] = cd
        airfoil_data[j, 3] = cm

    polar_data_test1 = ["polar_data", airfoil_data]
    polar_data_test2 = ["polar_data", airfoil_data * 1.1]  # Slightly different data

    # Create two sections with slightly different polar data
    section1 = MockSection([0, 0, 0], [1, 0, 0], polar_data_test1)
    section2 = MockSection([0, 10, 0], [1, 10, 0], polar_data_test2)

    # Create panel
    panel = create_panel(
        section1,
        section2,
    )

    assert hasattr(panel, "_panel_polar_data")
    assert panel._panel_polar_data is not None
    assert panel._panel_polar_data.shape == airfoil_data.shape

    # Check if panel_polar_data is correctly averaged
    expected_data = (airfoil_data + airfoil_data * 1.1) / 2
    assert np.allclose(panel._panel_polar_data, expected_data, atol=1e-6)


def test_panel_aerodynamic_center(sample_panel):
    expected_ac = np.array([0.25, 5, 0])
    assert np.allclose(sample_panel.aerodynamic_center, expected_ac)


def test_panel_control_point(sample_panel):
    expected_control_point = np.array([0.75, 5, 0])
    assert np.allclose(sample_panel.control_point, expected_control_point)


def test_panel_reference_frame(sample_panel):
    # Calculate the local reference frame
    # x_airf defined upwards from the chord-line, perpendicular to the panel
    # y_airf defined parallel to the chord-line, from LE-to-TE
    # z_airf along the LE, in plane (towards left tip, along span) from the airfoil perspective

    # you can think about subtracting of vectors like:
    # to get from A to B you do C = B - A, so work from backwards.

    LE_point_1 = sample_panel.LE_point_1
    LE_point_2 = sample_panel.LE_point_2
    TE_point_1 = sample_panel.TE_point_1
    TE_point_2 = sample_panel.TE_point_2

    print(f"LE_point_1: {LE_point_1}")
    print(f"LE_point_2: {LE_point_2}")
    print(f"TE_point_1: {TE_point_1}")
    print(f"TE_point_2: {TE_point_2}")
    print(f"rib_1 (LE - TE): {LE_point_1 - TE_point_1}")
    print(f"rib_2 (LE - TE): {LE_point_2 - TE_point_2}")
    print(f"rib_1 (TE - LE): {TE_point_1 - LE_point_1}")
    print(f"rib_2 (TE - LE): {TE_point_2 - LE_point_2}")

    mid_LE_point = LE_point_1 + 0.5 * (LE_point_2 - LE_point_1)
    mid_TE_point = TE_point_1 + 0.5 * (TE_point_2 - TE_point_1)
    vec_LE_to_TE = mid_TE_point - mid_LE_point

    print(f"mid_LE_point: {mid_LE_point}")
    print(f"mid_TE_point: {mid_TE_point}")
    print(f"vec_LE_to_TE: {vec_LE_to_TE}")

    y_airf = vec_LE_to_TE / np.linalg.norm(vec_LE_to_TE)
    bound_point_1 = LE_point_1 + 0.25 * (TE_point_1 - LE_point_1)
    bound_point_2 = LE_point_2 + 0.25 * (TE_point_2 - LE_point_2)

    print(f"section_1_aerodynamic_center/bound_point_1: {bound_point_1}")
    print(f"section_2_aerodynamic_center/bound_point_2: {bound_point_2}")
    z_airf = (bound_point_2 - bound_point_1) / np.linalg.norm(
        bound_point_2 - bound_point_1
    )
    x_airf = np.cross(y_airf, z_airf)

    print(f"vec_LE_to_TE: {vec_LE_to_TE}")
    print(f"x_airf: {x_airf}")
    print(f"y_airf: {y_airf}")
    print(f"z_airf: {z_airf}")

    # testing against expected values
    assert np.allclose(sample_panel.x_airf, [0, 0, 1])
    assert np.allclose(sample_panel.y_airf, [1, 0, 0])
    assert np.allclose(sample_panel.z_airf, [0, 1, 0])

    # testing against algorithm
    assert np.allclose(sample_panel.x_airf, x_airf)
    assert np.allclose(sample_panel.y_airf, y_airf)
    assert np.allclose(sample_panel.z_airf, z_airf)

    # testing the bound_points
    assert np.allclose(sample_panel.bound_point_1, bound_point_1)
    assert np.allclose(sample_panel.bound_point_2, bound_point_2)


def test_panel_custom_initialization():
    section1 = MockSection([1, 1, 1], [2, 1, 1], ["inviscid"])
    section2 = MockSection([1, 2, 1], [2, 2, 1], ["inviscid"])
    custom_panel = create_panel(section1, section2)

    assert np.allclose(custom_panel.aerodynamic_center, [1.25, 1.5, 1])
    assert np.allclose(custom_panel.control_point, [1.75, 1.5, 1])


def test_va_setter(sample_panel):
    test_va = np.array([1, 2, 3])
    sample_panel.va = test_va
    assert np.array_equal(sample_panel.va, test_va)


def test_calculate_relative_alpha_and_relative_velocity(sample_panel):
    sample_panel.va = np.array([10, 0, 0])
    induced_velocity = np.array([1, 1, 1])

    alpha_calc, relative_velocity_calc = (
        sample_panel.calculate_relative_alpha_and_relative_velocity(induced_velocity)
    )

    # Calculate terms of induced corresponding to the airfoil directions
    norm_airf = sample_panel.x_airf
    tan_airf = sample_panel.y_airf

    # Calculate relative velocity and angle of attack
    relative_velocity = sample_panel.va + induced_velocity
    vn = np.dot(norm_airf, relative_velocity)
    vtan = np.dot(tan_airf, relative_velocity)
    alpha = np.arctan(vn / vtan)

    print(f"vn: {vn}")
    print(f"vtan: {vtan}")
    print(f"alpha_calc: {alpha_calc}")
    print(f"alpha: {alpha}")
    print(f"relative_velocity_calc: {relative_velocity_calc}")
    print(f"relative_velocity: {relative_velocity}")
    assert np.isclose(alpha, alpha_calc)
    assert np.allclose(relative_velocity, relative_velocity_calc)


def test_calculate_velocity_induced_bound_2D(sample_panel):
    control_point = np.array([0.5, 5, 0])
    gamma = 1.0
    induced_velocity = sample_panel.calculate_velocity_induced_bound_2D(control_point)

    assert isinstance(induced_velocity, np.ndarray)
    assert induced_velocity.shape == (3,)  # 2D velocity


def test_velocity_induced_single_ring_semiinfinite(sample_panel):
    control_point = np.array([0.5, 5, 0])
    gamma = 1.0
    va_norm = 1
    va_unit = np.array([1, 0, 0])
    induced_velocity = sample_panel.calculate_velocity_induced_single_ring_semiinfinite(
        control_point, False, va_norm, va_unit, gamma, core_radius_fraction=0.01
    )

    assert isinstance(induced_velocity, np.ndarray)
    assert induced_velocity.shape == (3,)  # 3D velocity


def test_calculate_filaments_for_plotting(sample_panel):
    filaments_for_plotting = sample_panel.calculate_filaments_for_plotting()
    for filament in filaments_for_plotting:
        assert filament[0].shape == (3,)
        assert filament[1].shape == (3,)
        assert isinstance(filament[2], str)


# %% Testing cl,cd,cm calculation


def test_calculate_cl_and_cd_cm():
    # Generate mock polar data, using inviscid standards
    aoa = np.arange(-180, 180, 1)
    airfoil_data = np.empty((len(aoa), 4))
    for j, alpha in enumerate(aoa):
        cl, cd, cm = 2 * np.pi * np.deg2rad(alpha), 0.05, 0.01
        airfoil_data[j, 0] = np.deg2rad(alpha)
        airfoil_data[j, 1] = cl
        airfoil_data[j, 2] = cd
        airfoil_data[j, 3] = cm

    polar_data_test1 = ["polar_data", airfoil_data]

    # Create two sections with slightly different polar data
    inviscid_section1 = MockSection([0, 0, 0], [1, 0, 0], ["inviscid"])
    inviscid_section2 = MockSection([0, 10, 0], [1, 10, 0], ["inviscid"])
    polar_data_section1 = MockSection([0, 0, 0], [1, 0, 0], polar_data_test1)
    polar_data_section2 = MockSection([0, 10, 0], [1, 10, 0], polar_data_test1)
    lei_airfoil_section1 = MockSection(
        [0, 0, 0], [3, 0, 0], ["lei_airfoil_breukels", [0.12, 0.8]]
    )
    lei_airfoil_section2 = MockSection(
        [0, 10, 0], [3, 10, 0], ["lei_airfoil_breukels", [0.15, 0.7]]
    )

    # Create panels
    inviscid_panel_instance = create_panel(inviscid_section1, inviscid_section2)
    polar_data_panel_instance = create_panel(polar_data_section1, polar_data_section2)
    lei_airfoil_panel_instance = create_panel(
        lei_airfoil_section1, lei_airfoil_section2
    )
    # calculating average t and k
    t_avg = (0.12 + 0.15) / 2
    k_avg = (0.8 + 0.7) / 2
    t = t_avg
    k = k_avg

    # testing several angles
    test_alphas = [-10, 0, 10, 40]
    for alpha in test_alphas:
        alpha_rad = np.deg2rad(alpha)

        # inviscid panel
        cl_inviscid = inviscid_panel_instance.calculate_cl(alpha_rad)
        expected_cl_inviscid = 2 * np.pi * alpha_rad
        assert np.isclose(cl_inviscid, expected_cl_inviscid)

        cd_cm_inviscid = inviscid_panel_instance.calculate_cd_cm(alpha_rad)
        expected_cm_cd_inviscid = [0.0, 0.0]
        assert np.isclose(cd_cm_inviscid[0], expected_cm_cd_inviscid[1])

        # polar data panel
        cl_polar_data = polar_data_panel_instance.calculate_cl(alpha_rad)
        expected_cl_polar_data = 2 * np.pi * alpha_rad
        assert np.isclose(cl_polar_data, expected_cl_polar_data)

        cd_cm_polar_data = polar_data_panel_instance.calculate_cd_cm(alpha_rad)
        expected_cm_cd_polar_data = [0.01, 0.05]
        assert np.isclose(cd_cm_polar_data[0], expected_cm_cd_polar_data[1])

        # LEI airfoil panel
        cl_lei_airfoil = lei_airfoil_panel_instance.calculate_cl(alpha_rad)
        expected_cl_lei_airfoil = thesis_functions.LEI_airf_coeff(t, k, alpha)[0]
        print(f" ")
        print(f"alpha: {alpha}")
        print(f"expected_cl={expected_cl_lei_airfoil},calculating with t={t}, k={k}")
        print(f"cl_lei_airfoil: {cl_lei_airfoil}")
        assert np.isclose(cl_lei_airfoil, expected_cl_lei_airfoil)

        cd_cm_lei_airfoil = lei_airfoil_panel_instance.calculate_cd_cm(alpha_rad)
        expected_cd_cm_lei_airfoil = thesis_functions.LEI_airf_coeff(t, k, alpha)[1:]
        print(
            f"LEI Airfoil Cd/Cm: alpha={alpha}, cd={cd_cm_lei_airfoil[0]}, expected_cd={expected_cd_cm_lei_airfoil[0]}, cm={cd_cm_lei_airfoil[1]}, expected_cm={expected_cd_cm_lei_airfoil[1]}"
        )
        assert np.isclose(cd_cm_lei_airfoil[0], expected_cd_cm_lei_airfoil[0])
        assert np.isclose(cd_cm_lei_airfoil[1], expected_cd_cm_lei_airfoil[1])


def test_lei_airfoil_breukels_polynomial_new_against_old():

    # Create two sections with LEI airfoil parameters
    t1, k1 = 0.12, 0.8
    t2, k2 = 0.15, 0.7
    section1 = MockSection([0, 0, 0], [1, 0, 0], ["lei_airfoil_breukels", [t1, k1]])
    section2 = MockSection([0, 10, 0], [1, 10, 0], ["lei_airfoil_breukels", [t2, k2]])

    # Create panel
    panel = create_panel(section1, section2)

    # Check if LEI airfoil coefficients are correctly calculated
    assert hasattr(panel, "_cl_coefficients")
    assert hasattr(panel, "_cd_coefficients")
    assert hasattr(panel, "_cm_coefficients")

    # Test a few angles to ensure coefficients are calculated correctly
    test_angles = [-10, 0, 10]
    for alpha in test_angles:
        t = (t1 + t2) / 2
        k = (k1 + k2) / 2
        cl, cd, cm = thesis_functions.LEI_airf_coeff(t, k, alpha)

        calculated_cl = np.polyval(panel._cl_coefficients, alpha)
        calculated_cd = np.polyval(panel._cd_coefficients, alpha)
        calculated_cm = np.polyval(panel._cm_coefficients, alpha)

        assert np.isclose(calculated_cl, cl, atol=1e-6)
        assert np.isclose(calculated_cd, cd, atol=1e-6)
        assert np.isclose(calculated_cm, cm, atol=1e-6)


def get_v3_case_params():

    wing_type = "LEI_kite"
    dist = "lin"
    N_split = 4
    aoas = np.arange(-4, 24, 2)
    Umag = 22
    # convergence criteria
    max_iterations = 1500
    allowed_error = 1e-5
    relaxation_factor = 0.03
    core_radius_fraction = 1e-20

    # Wing geometry
    coord_struc = thesis_functions.get_CAD_matching_uri()
    coord = thesis_functions.struct2aero_geometry(coord_struc) / 1000

    N = len(coord) // 2

    # LE thickness at each section [m]
    # 10 sections
    LE_thicc = 0.1

    # Camber for each section (ct in my case)
    camber = 0.095

    # Refine structrural mesh into more panels
    coord = thesis_functions.refine_LEI_mesh(coord, N - 1, N_split)
    N = int(len(coord) / 2)  # Number of section after refining the mesh

    # Definition of airfoil coefficients
    # Based on Breukels (2011) correlation model
    aoas_for_polar = np.arange(-100, 100, 0.1)
    data_airf = np.empty((len(aoas_for_polar), 4))
    for j in range(len(aoas_for_polar)):
        alpha = aoas_for_polar[j]
        Cl, Cd, Cm = thesis_functions.LEI_airf_coeff(LE_thicc, camber, alpha)
        data_airf[j, 0] = alpha
        data_airf[j, 1] = Cl
        data_airf[j, 2] = Cd
        data_airf[j, 3] = Cm

    Atot = test_utils.calculate_projected_area(coord)
    coord_input_params = [coord, LE_thicc, camber]
    case_parameters = [
        coord_input_params,
        aoas,
        wing_type,
        Umag,
        0,
        Atot,
        max_iterations,
        allowed_error,
        relaxation_factor,
        core_radius_fraction,
        data_airf,
    ]

    return case_parameters


def test_lei_airfoil_breukels_polynomial_against_polar():

    (
        coord_input_params,
        aoas,
        wing_type,
        Umag,
        AR,
        Atot,
        max_iterations,
        allowed_error,
        relaxation_factor,
        core_radius_fraction,
        data_airf,
    ) = get_v3_case_params()
    [coord, LE_thicc, camber] = coord_input_params
    t = LE_thicc
    k = camber

    ### Creating polynomial and polar v3 panels
    # Polynomial
    aero_input_polynomial = ["lei_airfoil_breukels", [t, k]]
    section1 = MockSection([0, 0, 0], [1, 0, 0], aero_input_polynomial)
    section2 = MockSection([0, 1, 0], [1, 10, 0], aero_input_polynomial)
    panel_v3_polynomial = create_panel(section1, section2)
    # Polar
    # creating a data_airf in rad
    data_airf_rad = np.copy(data_airf)
    data_airf_rad[:, 0] = np.deg2rad(data_airf_rad[:, 0])
    aero_input_polar = ["polar_data", data_airf_rad]
    section1 = MockSection([0, 0, 0], [1, 0, 0], aero_input_polar)
    section2 = MockSection([0, 1, 0], [1, 10, 0], aero_input_polar)
    panel_v3_polar = create_panel(section1, section2)

    # TODO: somehow this only works when staying within the alpha range
    # of the interp data_airf, which maybe is not crazy
    # should check with old code if the data_airf was interpolated over a larger range
    # if so that should create bette results, closer to the same I hope

    test_angles = np.linspace(-50, 50, 5)
    for alpha in test_angles:
        alpha_rad = np.deg2rad(alpha)
        logging.info(f"--- alpha: {alpha}")
        ### Polynomial Calculation
        # old
        cl_polynomial_old, cd_polynomial_old, cm_polynomial_old = (
            thesis_functions.LEI_airf_coeff(t, k, alpha)
        )
        # new from panel function
        cl_polynomial_new = panel_v3_polynomial.calculate_cl(alpha_rad)
        cd_cm_polynomial_new = panel_v3_polynomial.calculate_cd_cm(alpha_rad)
        # asserting
        assert np.isclose(cl_polynomial_old, cl_polynomial_new)
        assert np.isclose(cd_polynomial_old, cd_cm_polynomial_new[0])
        assert np.isclose(cm_polynomial_old, cd_cm_polynomial_new[1])
        logging.info(f"cl_polynomial_old: {cl_polynomial_old}")
        logging.info(f"cl_polynomial_new: {cl_polynomial_new}")

        ### Polar Calculation
        # old
        cl_polar_old = np.interp(alpha, data_airf[:, 0], data_airf[:, 1])
        cd_polar_old = np.interp(alpha, data_airf[:, 0], data_airf[:, 2])
        cm_polar_old = np.interp(alpha, data_airf[:, 0], data_airf[:, 3])
        # new from panel function
        cl_polar_new = panel_v3_polar.calculate_cl(alpha_rad)
        cd_cm_polar_new = panel_v3_polar.calculate_cd_cm(alpha_rad)
        # asserting
        assert np.isclose(cl_polar_old, cl_polar_new)
        assert np.isclose(cd_polar_old, cd_cm_polar_new[0])
        assert np.isclose(cm_polar_old, cd_cm_polar_new[1])
        logging.info(f"cl_polar_old: {cl_polar_old}")
        logging.info(f"cl_polar_new: {cl_polar_new}")

        ### Comparing Polar and Polynomial
        assert np.isclose(cl_polar_new, cl_polynomial_new, atol=1e-7)
        assert np.isclose(cd_cm_polar_new[0], cd_cm_polynomial_new[0], atol=1e-7)
        assert np.isclose(cd_cm_polar_new[1], cd_cm_polynomial_new[1], atol=1e-7)
