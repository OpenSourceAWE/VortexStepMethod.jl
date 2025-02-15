# -*- coding: utf-8 -*-
import numpy as np
import logging
import sys
import os
import matplotlib.pyplot as plt
import pytest

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


def get_curved_case_params():

    max_chord = 2.18
    span = 6.969
    dist = "lin"
    N = 60
    aoas = np.arange(-4, 24, 1) / 180 * np.pi
    wing_type = "curved"
    Umag = 20
    AR = span / max_chord
    R = 4.673
    theta = 45 * np.pi / 180

    coord_input_params = [max_chord, span, theta, R, N, dist]
    # convergence criteria
    max_iterations = 1500
    allowed_error = 1e-5
    relaxation_factor = 0.03
    core_radius_fraction = 1e-20

    data_airf = np.loadtxt(
        r"./polars/clarky_maneia.csv",
        delimiter=",",
    )
    Atot = 14.40679
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


def test_curved():
    case_params = get_curved_case_params()
    # making sure not too many points are tested
    case_params[1] = np.deg2rad(np.array([9]))
    # comparison solution
    aoas = case_params[1]
    AR = case_params[4]
    polars_Maneia = np.loadtxt(
        "./polars/curved_wing_polars_maneia.csv",
        delimiter=",",
    )
    alpha_Maneia = polars_Maneia[:, 0]
    CL_Maneia = polars_Maneia[:, 1]
    CD_Maneia = polars_Maneia[:, 2]
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
    CL_maneia_at_new_alphas = []
    CD_maneia_at_new_alphas = []
    for alpha in aoas:
        alpha = np.rad2deg(alpha)
        CL_maneia_at_new_alphas.append(np.interp(alpha, alpha_Maneia, CL_Maneia))
        CD_maneia_at_new_alphas.append(np.interp(alpha, alpha_Maneia, CD_Maneia))

    # checking new VSM close to Maneia
    assert np.allclose(CL_maneia_at_new_alphas, CL_VSM_new, atol=1e-1)
    assert np.allclose(CD_maneia_at_new_alphas, CD_VSM_new, atol=1e-2)


if __name__ == "__main__":

    case_params = get_curved_case_params()

    aoas = case_params[1]
    aoas = np.deg2rad(np.linspace(0, 20, 10))
    case_params[1] = aoas
    AR = case_params[4]
    # comparing solution
    polars_Maneia = np.loadtxt(
        "./polars/curved_wing_polars_maneia.csv",
        delimiter=",",
    )
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
        wing_type="curved",
        aoas=[polars_Maneia[:, 0], aoas],
        CL_list=[polars_Maneia[:, 1], CL_LLT, CL_LLT_new, CL_VSM, CL_VSM_new],
        CD_list=[polars_Maneia[:, 2], CD_LLT, CD_LLT_new, CD_VSM, CD_VSM_new],
        gamma_list=[gamma_LLT, gamma_LLT_new, gamma_VSM, gamma_VSM_new],
        labels=["RANS [Maneia]", "LLT", "LLT_new", "VSM", "VSM_new"],
    )
    labels = ["Polars Maneia", "LLT", "LLT_new", "VSM", "VSM_new"]
    CL_list = [polars_Maneia[:, 1], CL_LLT, CL_LLT_new, CL_VSM, CL_VSM_new]
    CD_list = [polars_Maneia[:, 2], CD_LLT, CD_LLT_new, CD_VSM, CD_VSM_new]
    for i, aoa in enumerate(aoas):
        print(f"aoa = {np.rad2deg(aoa)}")
        for label, CD, CL in zip(labels, CD_list, CL_list):
            print(f"{label}: CL = {CL}, CD = {CD}")

    plt.show()
