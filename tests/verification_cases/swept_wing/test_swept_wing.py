# -*- coding: utf-8 -*-
import numpy as np
import logging
import sys
import os
import matplotlib.pyplot as plt
import pytest
import pandas as pd

# Go back to root folder
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, root_path)
import tests.utils as test_utils


@pytest.fixture(autouse=True)
def change_test_dir(request):
    # Get the directory of the current test file
    test_dir = os.path.dirname(request.fspath)
    # Change the working directory to the test file's directory
    os.chdir(test_dir)


def get_swept_wing_case_params():

    N = 60
    #   Rectangular wing coordinates
    chord = np.ones(N) * 1
    gamma = 0 / 180 * np.pi
    AR = 12
    span = AR * max(chord)
    Atot = span * max(chord)
    beta = np.zeros(N)
    offset = np.tan(gamma) * span / 2
    # beta = np.concatenate((np.linspace(-np.pi/6,0,round(N/2)),np.linspace(0,-np.pi/6,round(N/2))),axis =None)
    twist = np.ones(N) * 0
    dist = "lin"

    gamma = np.arctan(offset / (span / 2)) * 180 / np.pi
    Atot = span * max(chord)

    Umag = 20

    polars_airf = pd.read_csv("./polars/NACA4415_XFOIL_Re3e6.csv", skiprows=9)
    alpha_airf = np.arange(-10, 30)
    data_airf = np.zeros((len(polars_airf), 4))
    data_airf[:, 0] = polars_airf["alpha"]
    data_airf[:, 1] = polars_airf["CL"]
    data_airf[:, 2] = polars_airf["CD"]

    polars_airf = np.loadtxt("./polars/NACA4415_CFD_Re3e6.csv", delimiter=",")
    data_airf = np.zeros((len(polars_airf), 4))
    data_airf[:, 0] = polars_airf[:, 0]
    data_airf[:, 1] = polars_airf[:, 1]

    coord_input_params = [chord, offset, span, twist, beta, N, dist]
    # aoas = np.arange(0, 45, 2) / 180 * np.pi
    aoas = np.arange(0, 20, 3) / 180 * np.pi
    wing_type = "swept_wing"
    max_iterations = 1500
    allowed_error = 1e-5
    relaxation_factor = 0.03
    core_radius_fraction = 1e-20

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


def test_swept_wing():
    case_params = get_swept_wing_case_params()
    # making sure not too many points are tested
    case_params[1] = np.deg2rad(np.array([3, 6, 9]))
    # comparison solution
    aoas = case_params[1]
    polars_CFD = np.loadtxt("./polars/0sweepAR12_CFD.csv", delimiter=",")
    alpha_CFD = polars_CFD[:, 0]
    CL_CFD = polars_CFD[:, 1]
    CD_CFD = np.zeros_like(CL_CFD)
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
    # checking LTT old close to LLT new
    assert np.allclose(CL_LLT, CL_LLT_new, atol=2e-2)
    assert np.allclose(CD_LLT, CD_LLT_new, atol=2e-3)

    # checking VSMs to be close to one another
    assert np.allclose(CL_VSM, CL_VSM_new, atol=4e-2)
    assert np.allclose(CD_VSM, CD_VSM_new, atol=4e-3)

    # checking the LLT to be close to the VSM, with HIGHER tolerance
    tol_llt_to_vsm_CL = 2e-1
    tol_llt_to_vsm_CD = 4e-2
    assert np.allclose(CL_LLT, CL_VSM, atol=tol_llt_to_vsm_CL)
    assert np.allclose(CD_LLT, CD_VSM, atol=tol_llt_to_vsm_CD)
    assert np.allclose(CL_LLT_new, CL_VSM_new, atol=tol_llt_to_vsm_CL)
    assert np.allclose(CD_LLT_new, CD_VSM_new, atol=tol_llt_to_vsm_CD)

    # interpolating Maneia results to match the alphas
    CL_CFD_at_new_alphas = []
    CD_CFD_at_new_alphas = []
    for alpha in aoas:
        alpha = np.rad2deg(alpha)
        CL_CFD_at_new_alphas.append(np.interp(alpha, alpha_CFD, CL_CFD))
        CD_CFD_at_new_alphas.append(np.interp(alpha, alpha_CFD, CD_CFD))

    # checking new VSM close to CFD
    assert np.allclose(CL_CFD_at_new_alphas, CL_VSM_new, atol=1e-1)


if __name__ == "__main__":

    case_params = get_swept_wing_case_params()

    aoas = case_params[1]
    AR = case_params[4]
    # comparing solution
    polars_CFD = np.loadtxt("./polars/0sweepAR12_CFD.csv", delimiter=",")
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

    test_utils.plotting_CL_CD_gamma_LLT_VSM_old_new_comparison(
        panel_y=panel_y,
        AR=AR,
        wing_type="swept_wing",
        aoas=[polars_CFD[:, 0], aoas],
        CL_list=[polars_CFD[:, 1], CL_LLT, CL_LLT_new, CL_VSM, CL_VSM_new],
        CD_list=[
            np.zeros_like(polars_CFD[:, 0]),
            CD_LLT,
            CD_LLT_new,
            CD_VSM,
            CD_VSM_new,
        ],
        gamma_list=[gamma_LLT, gamma_LLT_new, gamma_VSM, gamma_VSM_new],
        labels=["RANS [Maneia]", "LLT", "LLT_new", "VSM", "VSM_new"],
    )
    labels = ["polars_CFD", "LLT", "LLT_new", "VSM", "VSM_new"]
    CL_list = [polars_CFD[:, 1], CL_LLT, CL_LLT_new, CL_VSM, CL_VSM_new]
    CD_list = [np.zeros_like(polars_CFD[:, 0]), CD_LLT, CD_LLT_new, CD_VSM, CD_VSM_new]

    for i, aoa in enumerate(aoas):
        print(f"aoa = {np.rad2deg(aoa)}")
        for label, CD, CL in zip(labels, CD_list, CL_list):
            print(f"{label}: CL = {CL}, CD = {CD}")

    plt.show()
