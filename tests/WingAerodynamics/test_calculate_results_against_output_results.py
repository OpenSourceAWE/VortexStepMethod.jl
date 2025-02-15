import pytest
import numpy as np
import logging
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
    flip_created_coord_in_pairs,
    calculate_projected_area,
)
import tests.thesis_functions_oriol_cayon as thesis_functions


def test_calculate_results():
    # Setup
    density = 1.225
    N = 40
    max_chord = 1
    span = 15.709  # AR = 20
    # span = 2.36  # AR = 3
    Umag = 20
    AR = span**2 / (np.pi * span * max_chord / 4)
    aoa = np.deg2rad(5)
    Uinf = np.array([np.cos(aoa), 0, np.sin(aoa)]) * Umag
    model = "VSM"

    # ### OLD FROM GEOMETRIC INPUT ####
    # dist = "cos"
    # coord = thesis_functions.generate_coordinates_el_wing(max_chord, span, N, dist)
    # Atot = max_chord / 2 * span / 2 * np.pi
    # ring_geo = "5fil"
    # alpha_airf = np.arange(-10, 30)
    # data_airf = np.zeros((len(alpha_airf), 4))
    # data_airf[:, 0] = alpha_airf
    # data_airf[:, 1] = alpha_airf / 180 * np.pi * 2 * np.pi
    # data_airf[:, 2] = alpha_airf * 0
    # data_airf[:, 3] = alpha_airf * 0
    # Gamma0 = np.zeros(N - 1)
    # conv_crit = {"Niterations": 1500, "error": 1e-5, "Relax_factor": 0.05}
    # # Define system of vorticity
    # controlpoints, rings, bladepanels, ringvec, coord_L = (
    #     thesis_functions.create_geometry_general(coord, Uinf, N, ring_geo, model)
    # )
    # # Solve for Gamma
    # Fmag, Gamma, aero_coeffs = (
    #     thesis_functions.solve_lifting_line_system_matrix_approach_semiinfinite(
    #         ringvec, controlpoints, rings, Uinf, Gamma0, data_airf, conv_crit, model
    #     )
    # )
    # # Calculate results using the reference function
    # F_rel_ref, F_gl_ref, Ltot_ref, Dtot_ref, CL_ref, CD_ref, CS_ref = (
    #     thesis_functions.output_results(
    #         Fmag, aero_coeffs, ringvec, Uinf, controlpoints, Atot
    #     )
    # )

    ### NEW ####
    dist = "cos"
    core_radius_fraction = 1e-20
    coord = thesis_functions.generate_coordinates_el_wing(max_chord, span, N, dist)
    coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
    wing = Wing(N, "unchanged")
    for idx in range(int(len(coord_left_to_right) / 2)):
        logging.debug(f"coord_left_to_right[idx] = {coord_left_to_right[idx]}")
        wing.add_section(
            coord_left_to_right[2 * idx],
            coord_left_to_right[2 * idx + 1],
            ["inviscid"],
        )
    wing_aero = WingAerodynamics([wing])
    wing_aero.va = Uinf

    ### Running the analysis
    solver_object = Solver(
        aerodynamic_model_type=model, core_radius_fraction=core_radius_fraction
    )
    # Solve the aerodynamics
    results_NEW = solver_object.solve(wing_aero)

    # Check the type and structure of the output
    assert isinstance(results_NEW, dict), "calculate_results should return a dictionary"

    ### OLD FROM WING OUTPUT ###
    # Calculating Fmag, using UNCORRECTED alpha
    alpha = results_NEW["alpha_uncorrected"]
    dyn_visc = 0.5 * density * np.linalg.norm(Uinf) ** 2
    n_panels = len(wing_aero.panels)
    lift, drag, moment = np.zeros(n_panels), np.zeros(n_panels), np.zeros(n_panels)
    for i, panel in enumerate(wing_aero.panels):
        lift[i] = dyn_visc * panel.calculate_cl(alpha[i]) * panel.chord
        drag[i] = dyn_visc * panel.calculate_cd_cm(alpha[i])[0] * panel.chord
        moment[i] = dyn_visc * panel.calculate_cd_cm(alpha[i])[1] * (panel.chord**2)
        print("lift:", lift, "drag:", drag, "moment:", moment)
    Fmag = np.column_stack([lift, drag, moment])

    # Calculating aero_coeffs, using CORRECTED alpha
    alpha = results_NEW["alpha_at_ac"]
    aero_coeffs = np.column_stack(
        (
            [alpha[i] for i, panel in enumerate(wing_aero.panels)],
            [panel.calculate_cl(alpha[i]) for i, panel in enumerate(wing_aero.panels)],
            [
                panel.calculate_cd_cm(alpha[i])[0]
                for i, panel in enumerate(wing_aero.panels)
            ],
            [
                panel.calculate_cd_cm(alpha[i])[1]
                for i, panel in enumerate(wing_aero.panels)
            ],
        )
    )
    ringvec = [{"r0": panel.width * panel.z_airf} for panel in wing_aero.panels]
    controlpoints = [
        {"tangential": panel.y_airf, "normal": panel.x_airf}
        for panel in wing_aero.panels
    ]
    Atot = calculate_projected_area(coord)

    # Calculate results using the reference function
    F_rel_ref, F_gl_ref, Ltot_ref, Dtot_ref, CL_ref, CD_ref, CS_ref = (
        thesis_functions.output_results(
            Fmag, aero_coeffs, ringvec, Uinf, controlpoints, Atot
        )
    )

    # Debug: Print the compared results
    cl_calculated = results_NEW["cl"]
    cd_calculated = results_NEW["cd"]
    cs_calculated = results_NEW["cs"]
    L_calculated = results_NEW["lift"]
    D_calculated = results_NEW["drag"]

    logging.info(f"cl_calculated: {cl_calculated}, CL_ref: {CL_ref}")
    logging.info(f"cd_calculated: {cd_calculated}, CD_ref: {CD_ref}")
    logging.info(f"cs_calculated: {cs_calculated}, CS_ref: {CS_ref}")
    logging.info(f"L_calculated: {L_calculated}, Ltot_ref: {Ltot_ref}")
    logging.info(f"D_calculated: {D_calculated}, Dtot_ref: {Dtot_ref}")

    ##########################
    ### COMPARING
    ##########################

    # Assert that the results are close
    np.testing.assert_allclose(cl_calculated, CL_ref, rtol=1e-4)
    np.testing.assert_allclose(cd_calculated, CD_ref, rtol=1e-4)
    np.testing.assert_allclose(cs_calculated, CS_ref, rtol=1e-4)
    np.testing.assert_allclose(L_calculated, Ltot_ref, rtol=1e-4)
    np.testing.assert_allclose(D_calculated, Dtot_ref, rtol=1e-4)

    # Check the shape of array outputs
    assert len(results_NEW["cl_distribution"]) == len(wing_aero.panels)
    assert len(results_NEW["cd_distribution"]) == len(wing_aero.panels)
