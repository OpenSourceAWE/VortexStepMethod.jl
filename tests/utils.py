import numpy as np
import logging
import time
import os
import sys
from copy import deepcopy
import matplotlib.pyplot as plt
import seaborn as sns

# Go back to root folder
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, root_path)
import tests.thesis_functions_oriol_cayon as thesis_functions

from VSM.WingGeometry import Wing
from VSM.WingAerodynamics import WingAerodynamics
from VSM.Solver import Solver


def cosspace(min, max, n_points):
    """
    Create an array with cosine spacing, from min to max values, with n points
    """
    mean = (max + min) / 2
    amp = (max - min) / 2
    return mean + amp * np.cos(np.linspace(np.pi, 0, n_points))


def generate_coordinates_rect_wing(chord, span, twist, beta, N, dist):
    """
    Generates the 3D coordinates of a rectangular wing with twist and dihedral.

    Parameters:
    chord: array-like
        The chord lengths of the wing panels.
    span: float
        The total span of the wing.
    twist: array-like
        The twist angles of the wing panels, in radians.
    beta: array-like
        The dihedral angles of the wing panels, in radians.
    N: int
        The number of panels along the span of the wing.
    dist: str
        The distribution type of the spanwise positions ('cos' for cosine-spaced, 'lin' for linear).

    Returns:
    coord: np.ndarray
        A 2N x 3 array of 3D coordinates for the leading and trailing edges of the wing panels.
    """
    coord = np.empty((2 * N, 3))
    if dist == "cos":
        span = cosspace(-span / 2, span / 2, N)
    elif dist == "lin":
        span = np.linspace(-span / 2, span / 2, N)

    for i in range(N):
        coord[2 * i, :] = np.array(
            [
                -0 * chord[i] * np.cos(twist[i]),
                span[i],
                0 * chord[i] * np.sin(twist[i]) - abs(span[i] * np.sin(beta[i])),
            ]
        )
        coord[2 * i + 1, :] = np.array(
            [
                1 * chord[i] * np.cos(twist[i]),
                span[i],
                -1 * chord[i] * np.sin(twist[i]) - abs(span[i] * np.sin(beta[i])),
            ]
        )

    return coord


def generate_coordinates_curved_wing(chord, span, theta, R, N, dist):
    """
    theta:
        This represents the angular extent of the curvature of the wing, in radians.
        It defines the range of angles over which the wing coordinates are distributed.
        A larger theta means a more pronounced curvature.

    R:
        This represents the radius of curvature of the wing.
        It defines the distance from the center of curvature (the origin in this case)
        to the points on the curved surface of the wing.
        A larger R means a larger radius, resulting in a more gradual curve.
    """
    coord = np.empty((2 * N, 3))
    if dist == "cos":
        theta = cosspace(-theta, theta, N)
    elif dist == "lin":
        theta = np.linspace(-theta, theta, N)
    elif dist == "cos2":
        theta1 = cosspace(-theta, -theta / N / 10, int(N / 2))
        theta2 = cosspace(theta / N / 10, theta, int(N / 2))
        theta = np.concatenate((theta1, theta2))

    for i in range(N):
        coord[2 * i, :] = np.array([0, R * np.sin(theta[i]), R * np.cos(theta[i])])
        coord[2 * i + 1, :] = np.array(
            [chord, R * np.sin(theta[i]), R * np.cos(theta[i])]
        )

    return coord


def generate_coordinates_el_wing(max_chord, span, N, dist):
    # cosine van garrel
    coord = np.empty((2 * N, 3))
    start = span * 1e-5
    if dist == "cos":
        y_arr = cosspace(-span / 2 + start, span / 2 - start, N)
    elif dist == "lin":
        y_arr = np.linspace(-span / 2 + start, span / 2 - start, N)

    c_arr = 2 * np.sqrt(1 - (y_arr / (span / 2)) ** 2) * max_chord / 2

    for i in range(N):
        coord[2 * i, :] = [-0.25 * c_arr[i], y_arr[i], 0]
        coord[2 * i + 1, :] = [0.75 * c_arr[i], y_arr[i], 0]

    return coord


# Check if the results are the same in list of dictionaries
def asserting_all_elements_in_list_dict(variable1, variable_expected):
    for i, (variable1_i, variable_expected_i) in enumerate(
        zip(variable1, variable_expected)
    ):
        for j, (key1, key2) in enumerate(
            zip(variable1_i.keys(), variable_expected_i.keys())
        ):
            logging.debug(f"key1 {key1}, key2 {key2}")
            logging.debug(f"variable1_i {variable1_i[key1]}")
            logging.debug(f"variable_expected_i {variable_expected_i[key1]}")
            # assert key1 == key2
            assert np.allclose(variable1_i[key1], variable_expected_i[key1], atol=1e-5)


# Check if the results are the same in list of list of dictionaries
def asserting_all_elements_in_list_list_dict(variable1, variable_expected):
    for i, (list1, list_expected) in enumerate(zip(variable1, variable_expected)):
        logging.info(f"list1 {list1}")
        logging.info(f"list_expected {list_expected}")
        for j, (dict1, dict_expected) in enumerate(zip(list1, list_expected)):
            logging.info(f"dict1 {dict1}")
            logging.info(f"dict_expected {dict_expected}")
            assert dict1.keys() == dict_expected.keys()
            for key1, key2 in zip(dict1.keys(), dict_expected.keys()):
                logging.info(f"key1 {key1}, key2 {key2}")
                logging.info(f"dict1[key1] {dict1[key1]}")
                logging.info(f"dict_expected[key1] {dict_expected[key1]}")
                assert key1 == key2
                # check breaks when entry is a string
                if isinstance(dict1[key1], str):
                    assert dict1[key1] == dict_expected[key1]
                else:
                    assert np.allclose(dict1[key1], dict_expected[key1], atol=1e-5)


def create_controlpoints_from_wing_object(wing, model):
    result = []

    for panel in wing.panels:
        if model == "VSM":
            cp = {
                "coordinates": panel.control_point,
                "chord": panel.chord,
                "normal": panel.x_airf,
                "tangential": panel.y_airf,
                "arf_coord": np.column_stack(
                    [panel.x_airf, panel.y_airf, panel.z_airf]
                ),
                "coordinates_aoa": panel.aerodynamic_center,
            }
        elif model == "LLT":
            cp = {
                "coordinates": panel.aerodynamic_center,
                "chord": panel.chord,
                "normal": panel.x_airf,
                "tangential": panel.y_airf,
                "airf_coord": np.column_stack(
                    [panel.x_airf, panel.y_airf, panel.z_airf]
                ),
            }
        else:
            raise ValueError(f"Model {model} not recognized")

        result.append(cp)
    return result


def create_ring_from_wing_object(wing, gamma_data=None):
    result = []
    va_norm = wing.va / np.linalg.norm(wing.va)
    filaments_list = [panel.filaments for panel in wing.panels]
    if gamma_data is None:
        gamma_data = [0 for _ in filaments_list]
    else:
        gamma_data = [gamma_data for _ in filaments_list]
    for filaments, gamma in zip(filaments_list, gamma_data):
        ring_filaments = []

        # bound starts
        filament_bound = filaments[0]
        filament_trailing_1 = filaments[1]
        filament_trailing_2 = filaments[2]
        filament_semi_1 = filaments[3]
        filament_semi_2 = filaments[4]

        # appending them in the correct order
        ring_filaments.append(
            {
                "id": "bound",
                "x1": filament_bound.x1,
                "x2": filament_bound.x2,
                "Gamma": gamma,
            }
        )
        ring_filaments.append(
            {
                "x1": filament_trailing_1.x1,
                "x2": filament_trailing_1.x2,
                "Gamma": gamma,
                "id": "trailing1",
            }
        )
        ring_filaments.append(
            {
                "dir": va_norm,
                "id": "trailing_inf1",
                "x1": filament_semi_1.x1,
                "Gamma": gamma,
            }
        )
        # be wary of incorrect naming convention here.
        ring_filaments.append(
            {
                "x2": filament_trailing_2.x2,
                "x1": filament_trailing_2.x1,
                "Gamma": gamma,
                "id": "trailing1",
            }
        )
        ring_filaments.append(
            {
                "x1": filament_semi_2.x1,
                "dir": va_norm,
                "id": "trailing_inf2",
                "Gamma": gamma,
            }
        )

        result.append(ring_filaments)

    return result


def create_wingpanels_from_wing_object(wing):
    wingpanels = []

    coordinates = np.zeros((2 * (len(wing.panels) + 1), 3))
    n_panels = len(wing.panels)
    for i in range(n_panels):
        coordinates[2 * i] = wing.panels[i].LE_point_1
        coordinates[2 * i + 1] = wing.panels[i].TE_point_1
        coordinates[2 * i + 2] = wing.panels[i].LE_point_2
        coordinates[2 * i + 3] = wing.panels[i].TE_point_2

    # Go through all wing panels
    for i, panel in enumerate(wing.panels):

        # Identify points defining the panel
        section = {
            "p1": coordinates[2 * i, :],
            "p2": coordinates[2 * i + 2, :],
            "p3": coordinates[2 * i + 3, :],
            "p4": coordinates[2 * i + 1, :],
        }
        wingpanels.append(section)
    return wingpanels


def create_ring_vec_from_wing_object(wing, model):
    result = []

    for panel in wing.panels:
        bound_1 = panel.bound_point_1
        bound_2 = panel.bound_point_2
        if model == "VSM":
            evaluation_point = panel.control_point
        elif model == "LLT":
            evaluation_point = panel.aerodynamic_center

        # Calculate the required vectors
        logging.debug(f"bound_1 {bound_1}")
        logging.debug(f"bound_2 {bound_2}")
        r0 = bound_2 - bound_1
        r1 = evaluation_point - bound_1
        r2 = evaluation_point - bound_2
        r3 = evaluation_point - (bound_2 + bound_1) / 2

        # Add the calculated vectors to the result
        result.append({"r0": r0, "r1": r1, "r2": r2, "r3": r3})

    return result


def create_coord_L_from_wing_object(wing):
    coord_L = []
    for panel in wing.panels:
        coord_L.append(panel.aerodynamic_center)
    return coord_L


def create_geometry_from_wing_object(wing, model):
    controlpoints = create_controlpoints_from_wing_object(wing, model)
    rings = create_ring_from_wing_object(wing)
    wingpanels = create_wingpanels_from_wing_object(wing)
    ringvec = create_ring_vec_from_wing_object(wing, model)
    coord_L = create_coord_L_from_wing_object(wing)

    return controlpoints, rings, wingpanels, ringvec, coord_L


def flip_created_coord_in_pairs(coord):
    # Reshape the array into pairs
    reshaped = coord.reshape(-1, 2, 3)

    # Reverse the order of the pairs
    flipped = np.flip(reshaped, axis=0)

    # Flatten back to the original shape
    return flipped.reshape(-1, 3)


def print_matrix(matrix, name="Matrix"):
    # Use np.array2string to convert the matrix to a nicely formatted string
    matrix_str = np.array2string(matrix, formatter={"float_kind": lambda x: "%.3f" % x})
    print(f"{name}:\n{matrix_str}")


def calculate_old_for_alpha_range(case_params):
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
    ) = case_params

    ## defining coord
    if wing_type == "elliptical":
        max_chord, span, N, dist = coord_input_params
        coord = thesis_functions.generate_coordinates_el_wing(max_chord, span, N, dist)
    elif wing_type == "curved":
        chord, span, theta, R, N, dist = coord_input_params
        coord = thesis_functions.generate_coordinates_curved_wing(
            chord, span, theta, R, N, dist
        )
    elif wing_type == "swept_wing":
        chord, offset, span, twist, beta, N, dist = coord_input_params
        coord = thesis_functions.generate_coordinates_swept_wing(
            chord, offset, span, twist, beta, N, dist
        )
    elif wing_type == "LEI_kite" or wing_type == "LEI_kite_polars":
        coord, t, k = coord_input_params
        N = len(coord) // 2
    else:
        raise ValueError(f"wing_type: {wing_type} not recognized")

    # aoa = 5.7106 * np.pi / 180
    # Uinf = np.array([np.cos(aoa), 0, np.sin(aoa)]) * Umag
    # Uinf = np.array([np.sqrt(0.99),0,0.1])
    conv_crit = {
        "Niterations": max_iterations,
        "error": allowed_error,
        "Relax_factor": relaxation_factor,
    }

    ring_geo = "5fil"
    model = "VSM"

    ## SOLVER + OUTPUT F
    start_time = time.time()
    CL1 = np.zeros(len(aoas))
    CL2 = np.zeros(len(aoas))
    CD1 = np.zeros(len(aoas))
    CD2 = np.zeros(len(aoas))
    Gamma_LLT = []
    Gamma_VSM = []
    for i in range(len(aoas)):

        Gamma0 = np.zeros(N - 1)

        Uinf = np.array([np.cos(aoas[i]), 0, np.sin(aoas[i])]) * Umag
        model = "LLT"
        # Define system of vorticity
        controlpoints, rings, bladepanels, ringvec, coord_L = (
            thesis_functions.create_geometry_general(coord, Uinf, N, ring_geo, model)
        )
        # Solve for Gamma
        Fmag, Gamma, aero_coeffs = (
            thesis_functions.solve_lifting_line_system_matrix_approach_semiinfinite(
                ringvec, controlpoints, rings, Uinf, Gamma0, data_airf, conv_crit, model
            )
        )
        # Output forces
        F_rel, F_gl, Ltot, Dtot, CL1[i], CD1[i], CS = thesis_functions.output_results(
            Fmag, aero_coeffs, ringvec, Uinf, controlpoints, Atot
        )
        Gamma_LLT.append(Gamma)
        model = "VSM"
        # Define system of vorticity
        controlpoints, rings, bladepanels, ringvec, coord_L = (
            thesis_functions.create_geometry_general(coord, Uinf, N, ring_geo, model)
        )
        # Solve for Gamma
        Fmag, Gamma, aero_coeffs = (
            thesis_functions.solve_lifting_line_system_matrix_approach_semiinfinite(
                ringvec, controlpoints, rings, Uinf, Gamma0, data_airf, conv_crit, model
            )
        )
        Gamma_VSM.append(Gamma)
        # Output forces
        F_rel, F_gl, Ltot, Dtot, CL2[i], CD2[i], CS = thesis_functions.output_results(
            Fmag, aero_coeffs, ringvec, Uinf, controlpoints, Atot
        )

        # Gamma0 = Gamma
        print(str((i + 1) / len(aoas) * 100) + " %")
    # end_time = time.time()
    # print(end_time - start_time)

    return CL1, CD1, CL2, CD2, Gamma_LLT, Gamma_VSM


def calculate_new_for_alpha_range(
    case_params,
    is_plotting=False,
):
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
    ) = case_params

    # transfering the data_airf first column to radians
    data_airf[:, 0] = np.deg2rad(data_airf[:, 0])

    ## defining coord
    if wing_type == "elliptical":
        max_chord, span, N, dist = coord_input_params
        coord = thesis_functions.generate_coordinates_el_wing(max_chord, span, N, dist)
        airfoil_input = ["inviscid"]
    elif wing_type == "curved":
        chord, span, theta, R, N, dist = coord_input_params
        coord = thesis_functions.generate_coordinates_curved_wing(
            chord, span, theta, R, N, dist
        )
        airfoil_input = ["polar_data", data_airf]
    elif wing_type == "swept_wing":
        chord, offset, span, twist, beta, N, dist = coord_input_params
        coord = thesis_functions.generate_coordinates_swept_wing(
            chord, offset, span, twist, beta, N, dist
        )
        airfoil_input = ["polar_data", data_airf]
    elif wing_type == "LEI_kite_polars":
        coord, thicc, camber = coord_input_params
        airfoil_input = ["polar_data", data_airf]
        N = len(coord) // 2
    elif wing_type == "LEI_kite":
        coord, thicc, camber = coord_input_params
        airfoil_input = ["lei_airfoil_breukels", [thicc, camber]]
        N = len(coord) // 2
    else:
        raise ValueError(f"wing_type: {wing_type} not recognized")

    coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
    wing = Wing(N, "unchanged")
    for idx in range(int(len(coord_left_to_right) / 2)):
        logging.debug(f"coord_left_to_right[idx] = {coord_left_to_right[idx]}")
        wing.add_section(
            coord_left_to_right[2 * idx],
            coord_left_to_right[2 * idx + 1],
            airfoil_input,
        )

    wing_aero = WingAerodynamics([wing])

    # initializing zero lists
    wing_aero = WingAerodynamics([wing])
    CL_LLT_new = np.zeros(len(aoas))
    CD_LLT_new = np.zeros(len(aoas))
    gamma_LLT_new = np.zeros((len(aoas), N - 1))
    CL_VSM_new = np.zeros(len(aoas))
    CD_VSM_new = np.zeros(len(aoas))
    gamma_VSM_new = np.zeros((len(aoas), N - 1))
    controlpoints_list = []
    LLT = Solver(
        aerodynamic_model_type="LLT",
        max_iterations=max_iterations,
        allowed_error=allowed_error,
        relaxation_factor=relaxation_factor,
        core_radius_fraction=core_radius_fraction,
    )
    # VSM
    VSM = Solver(
        aerodynamic_model_type="VSM",
        max_iterations=max_iterations,
        allowed_error=allowed_error,
        relaxation_factor=relaxation_factor,
        core_radius_fraction=core_radius_fraction,
    )
    for i, aoa_i in enumerate(aoas):
        logging.debug(f"aoa_i = {np.rad2deg(aoa_i)}")
        Uinf = np.array([np.cos(aoa_i), 0, np.sin(aoa_i)]) * Umag
        wing_aero.va = Uinf
        if i == 0 and is_plotting:
            raise ValueError("Plotting is not implemented for the first iteration")
        # LLT

        results_LLT = LLT.solve(wing_aero)
        CL_LLT_new[i] = results_LLT["cl"]
        CD_LLT_new[i] = results_LLT["cd"]
        gamma_LLT_new[i] = results_LLT["gamma_distribution"]

        results_VSM = VSM.solve(wing_aero)
        CL_VSM_new[i] = results_VSM["cl"]
        CD_VSM_new[i] = results_VSM["cd"]
        gamma_VSM_new[i] = results_VSM["gamma_distribution"]

        logging.debug(f"CD_LLT_new = {results_LLT['cd']}")
        logging.debug(f"CD_VSM_new = {results_VSM['cd']}")

        controlpoints_list.append(
            [panel.aerodynamic_center for panel in wing_aero.panels]
        )
    panel_y = [panel.aerodynamic_center[1] for panel in wing_aero.panels]
    AR = results_VSM["aspect_ratio_projected"]

    return (
        CL_LLT_new,
        CD_LLT_new,
        CL_VSM_new,
        CD_VSM_new,
        gamma_LLT_new,
        gamma_VSM_new,
        panel_y,
        AR,
    )


def plotting(
    x_axis_list: list,
    y_axis_list: list,
    labels: list,
    x_label: str,
    y_label: str,
    title: str,
    markers=None,
    alphas=None,
    colors=None,
    file_type=".pdf",
):

    if markers is None:
        markers = ["x", "x", ".", "."]
    if alphas is None:
        alphas = [0.8, 0.8, 0.8, 0.8]
    if colors is None:
        colors = sns.color_palette()

    plt_path = "./plots/"
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman"],
        }
    )
    plt.rcParams.update({"font.size": 10})
    plt.figure(figsize=(6, 4))
    for x_axis, y_axis, label, marker, alpha in zip(
        x_axis_list, y_axis_list, labels, markers, alphas
    ):
        plt.plot(x_axis, y_axis, marker=marker, alpha=alpha, label=label)
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.grid()
    plt.savefig(plt_path + title + file_type, bbox_inches="tight")


def plotting_CL_CD_gamma_LLT_VSM_old_new_comparison(
    panel_y, AR, wing_type, aoas, CL_list, CD_list, gamma_list, labels
):

    aoas_comparison = aoas[0]
    aoas = aoas[1]
    aoas_deg = np.rad2deg(aoas)
    plotting(
        x_axis_list=[aoas_comparison, aoas_deg, aoas_deg],
        y_axis_list=[CL_list[0], CL_list[1], CL_list[2]],
        labels=[labels[0], labels[1], labels[2]],
        x_label=r"$\alpha$ ($^\circ$)",
        y_label=r"$C_L$ ()",
        title=f"CL_alpha_LTT_{wing_type}_AR_" + str(round(AR, 1)),
        markers=None,
        alphas=None,
        colors=None,
        file_type=".pdf",
    )
    plotting(
        x_axis_list=[aoas_comparison, aoas_deg, aoas_deg],
        y_axis_list=[CD_list[0], CD_list[1], CD_list[2]],
        labels=[labels[0], labels[1], labels[2]],
        x_label=r"$\alpha$ ($^\circ$)",
        y_label=r"$C_D$ ()",
        title=f"CD_alpha_LTT_{wing_type}_AR_" + str(round(AR, 1)),
        markers=None,
        alphas=None,
        colors=None,
        file_type=".pdf",
    )
    plotting(
        x_axis_list=[aoas_comparison, aoas_deg, aoas_deg],
        y_axis_list=[CL_list[0], CL_list[3], CL_list[4]],
        labels=[labels[0], labels[3], labels[4]],
        x_label=r"$\alpha$ ($^\circ$)",
        y_label=r"$C_L$ ()",
        title=f"CL_alpha_VSM_{wing_type}_AR_" + str(round(AR, 1)),
        markers=None,
        alphas=None,
        colors=None,
        file_type=".pdf",
    )
    plotting(
        x_axis_list=[aoas_comparison, aoas_deg, aoas_deg],
        y_axis_list=[CD_list[0], CD_list[3], CD_list[4]],
        labels=[labels[0], labels[3], labels[4]],
        x_label=r"$\alpha$ ($^\circ$)",
        y_label=r"$C_D$ ()",
        title=f"CD_alpha_VSM_{wing_type}_AR_" + str(round(AR, 1)),
        markers=None,
        alphas=None,
        colors=None,
        file_type=".pdf",
    )
    # Plotting gamma for the mid aoa range
    idx = idx = int(len(aoas_deg) // 2)
    plotting(
        x_axis_list=[panel_y, panel_y, panel_y, panel_y],
        y_axis_list=[
            gamma_list[0][idx],
            gamma_list[1][idx],
            gamma_list[2][idx],
            gamma_list[3][idx],
        ],
        labels=[labels[1], labels[2], labels[3], labels[4]],
        x_label=r"$y$",
        y_label=r"$Gamma$",
        title=f"gamma_distribution_{wing_type}_AR_"
        + str(round(AR, 1))
        + "_aoa_"
        + str(round(AR, 1))
        + "_aoa_"
        + str(aoas_deg[idx]),
        markers=None,
        alphas=None,
        colors=None,
        file_type=".pdf",
    )


def calculate_projected_area(coord, z_plane_vector=np.array([0, 0, 1])):
    """Calculates the projected area of the wing onto a specified plane.

    The projected area is calculated based on the leading and trailing edge points of each section
    projected onto a plane defined by a normal vector (default is z-plane).

    Args:
        z_plane_vector (np.ndarray): Normal vector defining the projection plane (default is [0, 0, 1]).

    Returns:
        projected_area (float): The projected area of the wing.
    """
    # Normalize the z_plane_vector
    z_plane_vector = z_plane_vector / np.linalg.norm(z_plane_vector)

    # Helper function to project a point onto the plane
    def project_onto_plane(point, normal):
        return point - np.dot(point, normal) * normal

    projected_area = 0.0
    for i in range(len(coord) // 2 - 1):
        # Get the points for the current and next section
        LE_current = coord[2 * i]
        TE_current = coord[2 * i + 1]
        LE_next = coord[2 * (i + 1)]
        TE_next = coord[2 * (i + 1) + 1]

        # Project the points onto the plane
        LE_current_proj = project_onto_plane(LE_current, z_plane_vector)
        TE_current_proj = project_onto_plane(TE_current, z_plane_vector)
        LE_next_proj = project_onto_plane(LE_next, z_plane_vector)
        TE_next_proj = project_onto_plane(TE_next, z_plane_vector)

        # Calculate the lengths of the projected edges
        chord_current_proj = np.linalg.norm(TE_current_proj - LE_current_proj)
        chord_next_proj = np.linalg.norm(TE_next_proj - LE_next_proj)

        # Calculate the spanwise distance between the projected sections
        spanwise_distance_proj = np.linalg.norm(LE_next_proj - LE_current_proj)

        # Calculate the projected area of the trapezoid formed by these points
        area = 0.5 * (chord_current_proj + chord_next_proj) * spanwise_distance_proj
        projected_area += area

    return projected_area
