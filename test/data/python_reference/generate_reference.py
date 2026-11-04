# SPDX-FileCopyrightText: 2026 OpenSourceAWE contributors
# SPDX-License-Identifier: MPL-2.0
"""Write reference.yaml: inputs and results of the Python VSM v2.3.0 (awegroup/Vortex-Step-Method
at a2a91cb), read by test/verification/test_python_reference.jl.

Usage: python generate_reference.py <path to the Python repository> <output yaml>
"""
import sys
from pathlib import Path

import numpy as np
import yaml
from VSM.core.BodyAerodynamics import BodyAerodynamics
from VSM.core.Solver import Solver
from VSM.core.WingGeometry import Wing

ALPHA_GRID = np.deg2rad(np.arange(-10.0, 31.0, 2.0))
INVISCID = np.column_stack(
    (ALPHA_GRID, 2 * np.pi * ALPHA_GRID, 0 * ALPHA_GRID, 0 * ALPHA_GRID)
)


def stalling_polar():
    cl = np.where(
        ALPHA_GRID < np.deg2rad(14.0),
        2 * np.pi * ALPHA_GRID,
        2 * np.pi * np.deg2rad(14.0) - 1.5 * (ALPHA_GRID - np.deg2rad(14.0)),
    )
    return np.column_stack(
        (ALPHA_GRID, cl, 0.01 + 0.6 * ALPHA_GRID**2, -0.08 - 0.1 * ALPHA_GRID)
    )


def straight_sections():
    return [
        ([0.0, y, 0.0], [1.0, y, 0.0], INVISCID) for y in (5.0, -5.0)
    ]


def swept_sections():
    """Tapered, swept, twisted wing with dihedral."""
    sections = []
    for y in np.linspace(4.0, -4.0, 9):
        fraction = abs(y) / 4.0
        le = np.array([1.2 * fraction, y, 0.4 * fraction])
        chord = 1.5 - 0.7 * fraction
        twist = np.deg2rad(4.0 * (1.0 - fraction))
        te = le + chord * np.array([np.cos(twist), 0.0, -np.sin(twist)])
        sections.append((le.tolist(), te.tolist(), INVISCID))
    return sections


def arc_sections():
    """Circular-arc wing with anhedral and a tapered chord."""
    sections = []
    for angle in np.deg2rad(np.linspace(55.0, -55.0, 13)):
        le = np.array([0.3 * abs(angle), 4.0 * np.sin(angle), 4.0 * np.cos(angle)])
        chord = 1.2 - 0.4 * abs(angle)
        sections.append((le.tolist(), (le + [chord, 0.0, 0.0]).tolist(), stalling_polar()))
    return sections


def v3_sections(python_repo):
    """The V3 kite of the Python repository, with its Breukels polars on ALPHA_GRID."""
    path = (
        Path(python_repo)
        / "data/TUDELFT_V3_KITE/CAD_derived_geometry/aero_geometry_CAD_breukels_regression.yaml"
    )
    body = BodyAerodynamics.instantiate(n_panels=18, file_path=path)
    sections = []
    for section in body.wings[0].sections:
        polar = section.polar_data
        resampled = np.column_stack(
            [ALPHA_GRID] + [np.interp(ALPHA_GRID, polar[:, 0], polar[:, k]) for k in (1, 2, 3)]
        )
        sections.append(
            (list(map(float, section.LE_point)), list(map(float, section.TE_point)), resampled)
        )
    return sections


def plain(value):
    return np.asarray(value, dtype=float).tolist()


def run_case(name, sections, n_panels, va, keep_matrices=False, **options):
    wing = Wing(n_panels=n_panels, spanwise_panel_distribution="uniform")
    for le, te, polar in sections:
        wing.add_section(np.array(le), np.array(te), polar)
    body = BodyAerodynamics([wing])
    rates = options.pop("body_rates", None)
    pivot = options.pop("pivot", [0.0, 0.0, 0.0])
    va_array = np.asarray(va, dtype=float)
    if va_array.ndim == 2:
        va_array = np.array([va_array[0] + (va_array[1] - va_array[0]) * i / (n_panels - 1)
                             for i in range(n_panels)])
    if rates is None:
        body.va = va_array
    else:
        type(body).va.fset(body, va=va_array, body_rates=rates, body_axis=np.eye(3),
                           reference_point=np.array(pivot))
    reference_point = np.array([0.3, 0.1, 1.0])
    solver = Solver(allowed_error=1e-10, max_iterations=400000,
                    reference_point=reference_point, **options)
    results = solver.solve(body)
    assert results["gamma_converged"], name
    panels = body.panels
    va_panels = np.array([p.va for p in panels])
    speeds = np.linalg.norm(va_panels, axis=1)
    units = va_panels / speeds[:, None]
    model = options.get("aerodynamic_model_type", "VSM")
    AIC = np.array(body.compute_AIC_matrices(model, solver.core_radius_fraction, speeds, units))
    AIC_qc = np.array(
        body.compute_quarter_chord_AIC_matrices(solver.core_radius_fraction, speeds, units)
    )
    return {
        "name": name,
        "n_panels": n_panels,
        "sections": [{"LE": list(le), "TE": list(te), "polar": plain(polar)}
                     for le, te, polar in sections],
        "va": plain(va),
        "body_rates": None if rates is None else plain(rates),
        "pivot": plain(pivot),
        "reference_point": plain(reference_point),
        "options": options,
        "relaxation_factor": float(solver.relaxation_factor_used),
        "chord": plain([p.chord for p in panels]),
        "width": plain([p.width for p in panels]),
        "normal": plain([p.x_airf for p in panels]),
        "chordwise": plain([p.y_airf for p in panels]),
        "spanwise": plain([p.z_airf for p in panels]),
        "panel_va": plain(va_panels),
        "AIC": plain(AIC) if keep_matrices else None,
        "AIC_quarter_chord": plain(AIC_qc) if keep_matrices else None,
        "gamma": plain(results["gamma_distribution"]),
        "alpha": np.ravel(results["alpha_uncorrected"]).tolist(),
        "alpha_at_ac": np.ravel(results["alpha_at_ac"]).tolist(),
        "force": [float(results[key]) for key in ("Fx", "Fy", "Fz")],
        "moment": [float(results[key]) for key in ("Mx", "My", "Mz")],
        "force_coefficients": [float(results[key]) for key in ("cl", "cd", "cs")],
        "moment_coefficients": [float(results[key]) for key in ("cmx", "cmy", "cmz")],
        "lift_drag_side": [float(results[key]) for key in ("lift", "drag", "side")],
        "F_distribution": plain(results["F_distribution"]),
        "M_distribution": plain(results["M_distribution"]),
        "cl_distribution": np.ravel(results["cl_distribution"]).tolist(),
        "cd_distribution": np.ravel(results["cd_distribution"]).tolist(),
        "q_ref": float(results["q_ref"]),
        "va_ref": plain(results["va_ref"]),
        "projected_area": float(results["projected_area"]),
        "area_all_panels": float(results["area_all_panels"]),
    }


def main(python_repo, output):
    alpha = np.deg2rad(8.0)
    va = (15.0 * np.array([np.cos(alpha), 0.0, np.sin(alpha)])).tolist()
    sideslip = (15.0 * np.array([0.97, 0.15, 0.19])).tolist()
    arc = arc_sections()
    cases = [
        run_case("straight", straight_sections(), 20, va),
        run_case("swept", swept_sections(), 16, va, keep_matrices=True),
        run_case("arc", arc, 12, va, keep_matrices=True),
        run_case("arc_sideslip", arc, 12, sideslip),
        run_case("arc_llt", arc, 12, va, aerodynamic_model_type="LLT"),
        run_case("arc_viscous", arc, 12, va, is_with_viscous_drag_correction=True),
        run_case("arc_rates", arc, 12, va, body_rates=[0.2, -0.3, 0.5], pivot=[0.1, 0.0, 2.0]),
        run_case("arc_distributed", arc, 12, [[14.0, 0.0, 2.0], [17.0, 1.0, 2.5]]),
        run_case("v3", v3_sections(python_repo), 18, va),
    ]
    defaults = Solver()
    reference = {
        "python_commit": "a2a91cb02c05cef4c0201dfe2a666eba8e941f0f",
        "defaults": {
            "max_iterations": defaults.max_iterations,
            "core_radius_fraction": defaults.core_radius_fraction,
            "is_aoa_corrected": defaults.is_aoa_corrected,
            "is_with_attached_trailed_vortex_force": defaults.is_with_attached_trailed_vortex_force,
        },
        "cases": cases,
    }
    with open(output, "w") as file:
        yaml.safe_dump(reference, file, default_flow_style=None, width=100000)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
