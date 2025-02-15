import pytest
import numpy as np
import logging
import pprint
from copy import deepcopy
from VSM.Solver import Solver
from VSM.Panel import Panel
from VSM.WingAerodynamics import WingAerodynamics
from VSM.WingGeometry import Wing


import os
import sys

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, root_path)

from tests.utils import (
    generate_coordinates_el_wing,
    generate_coordinates_rect_wing,
    generate_coordinates_curved_wing,
    asserting_all_elements_in_list_dict,
    asserting_all_elements_in_list_list_dict,
    create_ring_from_wing_object,
    print_matrix,
    flip_created_coord_in_pairs,
)
from tests.thesis_functions_oriol_cayon import (
    vec_norm,
    dot_product,
    vector_projection,
    update_Gamma_single_ring,
    velocity_induced_single_ring_semiinfinite,
    velocity_induced_bound_2D,
    velocity_induced_single_ring_semiinfinite_nocore,
    create_geometry_general,
    update_Gamma_single_ring,
    velocity_3D_trailing_vortex_semiinfinite,
    velocity_3D_bound_vortex,
    velocity_3D_trailing_vortex,
)


def thesis_induction_matrix_creation(
    ringvec,
    controlpoints,
    rings,
    Uinf,
    Gamma0,
    data_airf,
    conv_crit,
    model,
):
    """
    Solve the VSM or LLM by finding the distribution of Gamma

    Parameters
    ----------
    ringvec : List of dictionaries containing the vectors that define each ring
    controlpoints :List of dictionaries with the variables needed to define each wing section
    rings : List of list with the definition of each vortex filament
    Uinf : Wind speed velocity vector
    data_airf : 2D airfoil data with alpha, Cl, Cd, Cm
    recalc_alpha : True if you want to recalculate the induced angle of attack at 1/4 of the chord (VSM)
    Gamma0 : Initial Guess of Gamma
    model : VSM: Vortex Step method/ LLT: Lifting Line Theory

    Returns
    -------
    F: Lift, Drag and Moment at each section
    Gamma: Gamma at each section
    aero_coeffs: alpha, cl, cd, cm at each wing section

    """

    nocore = False  # To shut down core corrections input True
    # Initialization of the parameters
    velocity_induced = []
    u = 0
    v = 0
    w = 0
    N = len(rings)
    Gammaini = Gamma0
    Gamma = np.zeros(N)
    GammaNew = Gammaini
    Lift = np.zeros(N)
    Drag = np.zeros(N)
    Ma = np.zeros(N)
    alpha = np.zeros(N)
    cl = np.zeros(N)
    cd = np.zeros(N)
    cm = np.zeros(N)
    MatrixU = np.empty((N, N))
    MatrixV = np.empty((N, N))
    MatrixW = np.empty((N, N))
    U_2D = np.zeros((N, 3))
    logging.debug(f"shape(U_2D): {np.shape(U_2D)}")

    # Number of iterations and convergence criteria
    Niterations = conv_crit["Niterations"]
    errorlimit = conv_crit["error"]
    ConvWeight = conv_crit["Relax_factor"]
    converged = False
    rho = 1.225

    coord_cp = [controlpoints[icp]["coordinates"] for icp in range(N)]
    chord = [controlpoints[icp]["chord"] for icp in range(N)]
    airf_coord = [controlpoints[icp]["airf_coord"] for icp in range(N)]

    for icp in range(N):

        if model == "VSM":
            # Velocity induced by a infinte bound vortex with Gamma = 1
            logging.debug(f"ringvec[icp]: {ringvec[icp]}")
            U_2D[icp] = velocity_induced_bound_2D(ringvec[icp])

        for jring in range(N):
            rings[jring] = update_Gamma_single_ring(rings[jring], 1, 1)
            # Calculate velocity induced by a ring to a control point
            velocity_induced = velocity_induced_single_ring_semiinfinite(
                rings[jring], coord_cp[icp], model, vec_norm(Uinf)
            )
            # If CORE corrections are deactivated
            if nocore == True:
                # Calculate velocity induced by a ring to a control point
                velocity_induced = velocity_induced_single_ring_semiinfinite_nocore(
                    rings[jring], coord_cp[icp], model
                )

            # AIC Matrix
            MatrixU[icp, jring] = velocity_induced[0]
            MatrixV[icp, jring] = velocity_induced[1]
            MatrixW[icp, jring] = velocity_induced[2]

            # THIS is DIFFERENT from the old thesis code
            # as this was considered wrong in that code
            if icp == jring:
                if model == "VSM":
                    logging.debug(f"U_2D[icp]: {U_2D[icp]}")
                    MatrixU[icp, jring] -= U_2D[icp, 0]
                    MatrixV[icp, jring] -= U_2D[icp, 1]
                    MatrixW[icp, jring] -= U_2D[icp, 2]

    return MatrixU, MatrixV, MatrixW


def test_induction_matrix_creation():

    n_panels = 3
    # N = number of SECTIONS
    N = n_panels + 1
    max_chord = 1
    span = 2.36
    AR = span**2 / (np.pi * span * max_chord / 4)
    dist = "cos"
    coord = generate_coordinates_el_wing(max_chord, span, N, dist)

    Atot = max_chord / 2 * span / 2 * np.pi
    logging.debug(f"N: {N}")
    logging.debug(f"len(coord): {len(coord)}")

    Umag = 20
    aoa = 5.7106 * np.pi / 180
    Uinf = np.array([np.cos(aoa), 0, np.sin(aoa)]) * Umag
    # Uinf = np.array([np.sqrt(0.99),0,0.1])

    conv_crit = {"Niterations": 1500, "error": 1e-5, "Relax_factor": 0.05}

    Gamma0 = np.zeros(N - 1)

    ring_geo = "5fil"

    alpha_airf = np.arange(-10, 30)
    data_airf = np.zeros((len(alpha_airf), 4))
    data_airf[:, 0] = alpha_airf
    data_airf[:, 1] = alpha_airf / 180 * np.pi * 2 * np.pi
    data_airf[:, 2] = alpha_airf * 0
    data_airf[:, 3] = alpha_airf * 0

    # Creating the wing object
    ## Elliptical Wing
    coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
    core_radius_fraction = 1e-20  # only value I could find
    wing = Wing(n_panels, "unchanged")
    for idx in range(int(len(coord_left_to_right) / 2)):
        wing.add_section(
            coord_left_to_right[2 * idx], coord_left_to_right[2 * idx + 1], ["inviscid"]
        )
    wing_aero = WingAerodynamics([wing])
    wing_aero.va = Uinf

    ###########
    ### LTT ###
    ###########
    model = "LLT"

    ### Geometry ###
    ### THESIS ###
    # Define system of vorticity
    controlpoints, rings, bladepanels, ringvec, coord_L = create_geometry_general(
        coord, Uinf, N, ring_geo, model
    )
    MatrixU, MatrixV, MatrixW = thesis_induction_matrix_creation(
        deepcopy(ringvec),
        deepcopy(controlpoints),
        deepcopy(rings),
        deepcopy(Uinf),
        deepcopy(Gamma0),
        deepcopy(data_airf),
        deepcopy(conv_crit),
        deepcopy(model),
    )

    ### NEW ###
    va_norm_array = np.array([vec_norm(Uinf)] * len(coord))
    va_unit_array = np.array([Uinf / vec_norm(Uinf)] * len(coord))
    AIC_x, AIC_y, AIC_z = wing_aero.calculate_AIC_matrices(
        model, core_radius_fraction, va_norm_array, va_unit_array
    )

    print_matrix(MatrixU, "MatrixU")
    print_matrix(AIC_x, "AIC_x")
    print_matrix(MatrixV, "MatrixV")
    print_matrix(-AIC_y, "AIC_y")
    print_matrix(MatrixW, "MatrixW")
    print_matrix(AIC_z, "AIC_z")

    # Check if the matrices are the same
    assert np.allclose(MatrixU, AIC_x, atol=1e-5)
    assert np.allclose(MatrixV, -AIC_y, atol=1e-5)
    assert np.allclose(MatrixW, AIC_z, atol=1e-5)

    ###########
    ### VSM ###
    ###########

    model = "VSM"

    ### THESIS ###
    controlpoints, rings, bladepanels, ringvec, coord_L = create_geometry_general(
        coord, Uinf, N, ring_geo, model
    )
    MatrixU, MatrixV, MatrixW = thesis_induction_matrix_creation(
        deepcopy(ringvec),
        deepcopy(controlpoints),
        deepcopy(rings),
        deepcopy(Uinf),
        deepcopy(Gamma0),
        deepcopy(data_airf),
        deepcopy(conv_crit),
        deepcopy(model),
    )

    ### NEW ###

    # wing_aero.va = Uinf
    AIC_x, AIC_y, AIC_z = wing_aero.calculate_AIC_matrices(
        model, core_radius_fraction, va_norm_array, va_unit_array
    )
    # Check if the matrices are the same
    assert np.allclose(MatrixU, AIC_x, atol=1e-8)
    assert np.allclose(MatrixV, -AIC_y, atol=1e-8)
    assert np.allclose(MatrixW, AIC_z, atol=1e-8)
