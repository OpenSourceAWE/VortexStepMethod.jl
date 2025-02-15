# -*- coding: utf-8 -*-
import numpy as np
import logging
import sys
import os
import matplotlib.pyplot as plt

# Go back to root folder
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, root_path)
import tests.utils as test_utils


def get_elliptical_case_params():
    # general geometry
    max_chord = 1
    span = 15.709  # AR = 20
    # span = 2.36  # AR = 3
    dist = "cos"
    N = 40
    aoas = np.arange(0, 20, 1) / 180 * np.pi
    wing_type = "elliptical"

    coord_input_params = [max_chord, span, N, dist]
    # wind
    Umag = 20
    AR = span**2 / (np.pi * span * max_chord / 4)
    Atot = max_chord / 2 * span / 2 * np.pi
    # convergence criteria
    max_iterations = 1500
    allowed_error = 1e-5
    relaxation_factor = 0.05
    core_radius_fraction = 1e-20

    # data_airf
    alpha_airf = np.arange(-10, 30)
    data_airf = np.zeros((len(alpha_airf), 4))
    data_airf[:, 0] = alpha_airf
    data_airf[:, 1] = alpha_airf / 180 * np.pi * 2 * np.pi
    data_airf[:, 2] = alpha_airf * 0
    data_airf[:, 3] = alpha_airf * 0

    case_parameters = [
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
    ]
    return case_parameters


def test_elliptical():
    case_params = get_elliptical_case_params()
    # making sure not too many points are tested
    case_params[1] = np.deg2rad(np.array([3, 9]))
    # analytical solution
    aoas = case_params[1]
    AR = case_params[4]
    CL_th = 2 * np.pi * aoas / (1 + 2 / AR)
    CDi_th = CL_th**2 / np.pi / AR
    # OLD numerical
    CL_LLT, CD_LLT, CL_VSM, CD_VSM, gamma_LLT, gamma_VSM = (
        test_utils.calculate_old_for_alpha_range(case_params)
    )
    # NEW numerical
    (
        CL_LLT_new,
        CD_LLT_new,
        CL_VSM_new,
        CD_VSM_new,
        gamma_LLT_new,
        gamma_VSM_new,
        panel_y,
        AR_projected,
    ) = test_utils.calculate_new_for_alpha_range(
        case_params,
        is_plotting=False,
    )

    # checking all LLTs to be close
    assert np.allclose(CL_th, CL_LLT, atol=1e-2)
    assert np.allclose(CDi_th, CD_LLT, atol=1e-3)
    assert np.allclose(CL_th, CL_LLT_new, atol=1e-2)
    assert np.allclose(CDi_th, CD_LLT_new, atol=1e-3)
    assert np.allclose(gamma_LLT, gamma_LLT_new, atol=1e-2)

    # checking VSMs to be close to one another
    assert np.allclose(CL_VSM, CL_VSM_new, atol=1e-2)
    assert np.allclose(CD_VSM, CD_VSM_new, atol=1e-3)

    # checking the LLT to be close to the VSM, with HIGHER tolerance
    tol_llt_to_vsm_CL = 1e-1
    tol_llt_to_vsm_CD = 1e-3
    assert np.allclose(CL_th, CL_VSM, atol=tol_llt_to_vsm_CL)
    assert np.allclose(CDi_th, CD_VSM, atol=tol_llt_to_vsm_CD)
    assert np.allclose(CL_th, CL_VSM_new, atol=tol_llt_to_vsm_CL)
    assert np.allclose(CDi_th, CD_VSM_new, atol=tol_llt_to_vsm_CD)
    assert np.allclose(CL_LLT, CL_VSM, atol=tol_llt_to_vsm_CL)
    assert np.allclose(CD_LLT, CD_VSM, atol=tol_llt_to_vsm_CD)
    assert np.allclose(CL_LLT_new, CL_VSM_new, atol=tol_llt_to_vsm_CL)
    assert np.allclose(CD_LLT_new, CD_VSM_new, atol=tol_llt_to_vsm_CD)


if __name__ == "__main__":

    case_params = get_elliptical_case_params()

    # analytical solution
    aoas = case_params[1]
    AR = case_params[4]
    CL_th = 2 * np.pi * aoas / (1 + 2 / AR)
    CDi_th = CL_th**2 / np.pi / AR
    # OLD numerical
    CL_LLT, CD_LLT, CL_VSM, CD_VSM, gamma_LLT, gamma_VSM = (
        test_utils.calculate_old_for_alpha_range(case_params)
    )
    # NEW numerical
    (
        CL_LLT_new,
        CD_LLT_new,
        CL_VSM_new,
        CD_VSM_new,
        gamma_LLT_new,
        gamma_VSM_new,
        panel_y,
        AR_projected,
    ) = test_utils.calculate_new_for_alpha_range(
        case_params,
        is_plotting=False,
    )

    CL_list = [CL_th, CL_LLT, CL_LLT_new, CL_VSM, CL_VSM_new]
    CD_list = [CDi_th, CD_LLT, CD_LLT_new, CD_VSM, CD_VSM_new]
    labels = ["Analytic LLT", "LLT", "LLT_new", "VSM", "VSM_new"]
    test_utils.plotting_CL_CD_gamma_LLT_VSM_old_new_comparison(
        panel_y=panel_y,
        AR=AR,
        wing_type="curved",
        aoas=[aoas, aoas],
        CL_list=CL_list,
        CD_list=CD_list,
        gamma_list=[gamma_LLT, gamma_LLT_new, gamma_VSM, gamma_VSM_new],
        labels=labels,
    )
    for i, aoa in enumerate(aoas):
        print(f"aoa = {np.rad2deg(aoa)}")
        for label, CD, CL in zip(labels, CD_list, CL_list):
            print(f"{label}: CL = {CL[i]}, CD = {CD[i]}")
