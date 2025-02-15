import numpy as np
import logging
from copy import deepcopy
from VSM.Filament import BoundFilament, SemiInfiniteFilament
from VSM.WingGeometry import Wing
from VSM.WingAerodynamics import WingAerodynamics

import os
import sys

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, root_path)
from tests.utils import (
    generate_coordinates_el_wing,
    generate_coordinates_rect_wing,
    generate_coordinates_curved_wing,
    asserting_all_elements_in_list_list_dict,
    asserting_all_elements_in_list_dict,
    create_controlpoints_from_wing_object,
    create_ring_from_wing_object,
    create_wingpanels_from_wing_object,
    create_ring_vec_from_wing_object,
    create_coord_L_from_wing_object,
    flip_created_coord_in_pairs,
)
from tests.thesis_functions_oriol_cayon import create_geometry_general


def test_multiple():
    for model in ["VSM", "LLT"]:
        logging.debug(f"model: {model}")
        for wing_type in ["rectangular", "curved", "elliptical"]:
            logging.debug(f"wing_type: {wing_type}")
            wing_aero, coord, Uinf, model = create_geometry(model, wing_type)
            creating_tests(wing_aero, coord, Uinf, model)


def creating_tests(wing_aero, coord, Uinf, model):

    # Generate geometry
    (
        expected_controlpoints,
        expected_rings,
        expected_bladepanels,
        expected_ringvec,
        expected_coord_L,
    ) = create_geometry_general(coord, Uinf, int(len(coord) / 2), "5fil", model)

    for i, _ in enumerate(wing_aero.panels):
        logging.debug(f"i: {i}")
        #### handling of control points dict
        index_reversed = -(i + 1)
        panel = wing_aero.panels[index_reversed]
        if model == "VSM":
            evaluation_point = panel.control_point
        elif model == "LLT":
            evaluation_point = panel.aerodynamic_center

        assert np.allclose(
            evaluation_point, expected_controlpoints[i]["coordinates"], atol=1e-4
        )
        assert np.allclose(panel.chord, expected_controlpoints[i]["chord"], atol=1e-4)
        assert np.allclose(panel.x_airf, expected_controlpoints[i]["normal"], atol=1e-4)
        assert np.allclose(
            panel.y_airf, expected_controlpoints[i]["tangential"], atol=1e-4
        )
        assert np.allclose(
            np.column_stack((panel.x_airf, panel.y_airf, panel.z_airf)),
            expected_controlpoints[i]["airf_coord"],
            atol=1e-4,
        )
        if model == "VSM":
            assert np.allclose(
                panel.aerodynamic_center,
                expected_controlpoints[i]["coordinates_aoa"],
                atol=1e-4,
            )

        #### handling of rings dict
        expected_ring_i = expected_rings[i]
        expected_ring_i_list = [
            expected_ring_i[0],
            expected_ring_i[1],
            expected_ring_i[2],
            expected_ring_i[3],
            expected_ring_i[4],
        ]

        filaments = panel.filaments
        # they are differently arranged, so must reshuffled them
        filament_list = [
            filaments[0],
            filaments[2],
            filaments[4],
            filaments[1],
            filaments[3],
        ]
        for k, _ in enumerate(filament_list):
            logging.debug(f"exp_x1 {expected_ring_i_list[k]['x1']}")
        for k, _ in enumerate(filament_list):
            logging.debug(f"x1 {filament_list[k].x1}")
        for j, _ in enumerate(filament_list):
            logging.debug(f"j: {j}")
            if j == 0:  # bound filaments
                assert np.allclose(
                    filament_list[j].x1, expected_ring_i_list[j]["x1"], atol=1e-4
                )
                assert np.allclose(
                    filament_list[j].x2, expected_ring_i_list[j]["x2"], atol=1e-4
                )
            elif j == (1 or 3):  # trailing filaments
                assert np.allclose(
                    filament_list[j].x1, expected_ring_i_list[j]["x1"], atol=1e-4
                )
                assert np.allclose(
                    filament_list[j].x2, expected_ring_i_list[j]["x2"], atol=1e-4
                )
            else:  # semi-infinite filaments
                assert np.allclose(
                    filament_list[j].x1, expected_ring_i_list[j]["x1"], atol=1e-4
                )

        ### handling of bladepanels dict
        exp_bladepanels = expected_bladepanels[i]
        # rewriting these p1,2,3,4 points to match the reversed order
        p1 = panel.LE_point_2
        p2 = panel.LE_point_1
        p3 = panel.TE_point_1
        p4 = panel.TE_point_2
        assert np.allclose(p1, exp_bladepanels["p1"], atol=1e-4)
        assert np.allclose(p2, exp_bladepanels["p2"], atol=1e-4)
        assert np.allclose(p3, exp_bladepanels["p3"], atol=1e-4)
        assert np.allclose(p4, exp_bladepanels["p4"], atol=1e-4)

        ### handling of ringvec dict
        exp_ringvec = expected_ringvec[i]
        logging.debug(f"exp_ringvec {exp_ringvec}")
        bound_1 = panel.bound_point_1
        bound_2 = panel.bound_point_2

        r3 = evaluation_point - (bound_1 + bound_2) / 2
        r0 = bound_1 - bound_2
        assert np.allclose(r0, exp_ringvec["r0"], atol=1e-4)
        assert np.allclose(r3, exp_ringvec["r3"], atol=1e-4)

        ### handling of coord_L dict
        exp_coord_L = expected_coord_L[i]
        assert np.allclose(panel.aerodynamic_center, exp_coord_L)


def create_geometry(model="VSM", wing_type="rectangular", plotting=False, N=40):
    max_chord = 1
    span = 17
    AR = span**2 / (np.pi * span * max_chord / 4)
    print(f"AR: {AR}")
    Umag = 20
    aoa = 5.7106 * np.pi / 180
    Uinf = np.array([np.cos(aoa), 0, np.sin(aoa)]) * Umag

    if wing_type == "rectangular":
        twist = np.linspace(-0.5, 0.5, N)
        beta = np.linspace(-2, 2, N)
        coord = generate_coordinates_rect_wing(
            max_chord * np.ones(N),
            span,
            twist,
            beta,
            N=N,
            dist="lin",
        )
    elif wing_type == "curved":
        coord = generate_coordinates_curved_wing(
            max_chord, span, np.pi / 4, 5, N, "cos"
        )
    elif wing_type == "elliptical":
        coord = generate_coordinates_el_wing(max_chord, span, N, "cos")
    else:
        raise ValueError("Invalid wing type")

    coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
    wing = Wing(N, "unchanged")
    for i in range(int(len(coord_left_to_right) / 2)):
        wing.add_section(
            coord_left_to_right[2 * i], coord_left_to_right[2 * i + 1], ["inviscid"]
        )
    wing_aero = WingAerodynamics([wing])
    wing_aero.va = Uinf
    if plotting:
        wing_aero.plot()
    return wing_aero, coord, Uinf, model


if __name__ == "__main__":
    n_sections = 20
    wing_aero, coord, Uinf, model = create_geometry(
        model="VSM", wing_type="elliptical", plotting=True, N=n_sections
    )
    wing_aero, coord, Uinf, model = create_geometry(
        model="VSM", wing_type="curved", plotting=True, N=n_sections
    )
    wing_aero, coord, Uinf, model = create_geometry(
        model="VSM", wing_type="rectangular", plotting=True, N=n_sections
    )
