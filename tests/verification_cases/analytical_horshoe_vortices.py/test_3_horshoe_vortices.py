import numpy as np
import logging
import matplotlib.pyplot as plt
import pytest
import pprint
from VSM.Filament import SemiInfiniteFilament, BoundFilament


def calculate_induced_vel_3_horsheshoes():
    # Gamma values for the filaments
    Uinf = np.array([1, 0, 0])
    va_norm = np.linalg.norm(Uinf)
    va_unit = Uinf / va_norm

    # The Test Case
    gammas = [2, 10, 5]
    # evaluation_point = np.array([82, -1.1, 0])
    evaluation_point = np.array([0, 0, 0])

    # Analytical solutions for the combined induced velocities
    analytical_solutions = {
        2: np.array([0, 0, 0.1061]),
        10: np.array([0, 0, -1.5915]),
        5: np.array([0, 0, 0.2653]),
    }

    left_gamma2 = np.array([0, -1, 0])
    right_gamma2 = np.array([0, -3, 0])
    horseshoe1_filaments = [
        BoundFilament(x1=right_gamma2, x2=left_gamma2),
        SemiInfiniteFilament(
            np.array(left_gamma2),
            va_unit,
            va_norm,
            filament_direction=1,
        ),
        SemiInfiniteFilament(
            np.array(right_gamma2),
            va_unit,
            va_norm,
            filament_direction=-1,
        ),
    ]
    left_gamma10 = np.array([0, 1, 0])
    right_gamma10 = np.array([0, -1, 0])
    horseshoe2_filaments = [
        BoundFilament(x1=right_gamma10, x2=left_gamma10),
        SemiInfiniteFilament(
            np.array(left_gamma10),
            va_unit,
            va_norm,
            filament_direction=1,
        ),
        SemiInfiniteFilament(
            np.array(right_gamma10),
            va_unit,
            va_norm,
            filament_direction=-1,
        ),
    ]
    left_gamma5 = np.array([0, 3, 0])
    right_gamma5 = np.array([0, 1, 0])
    horseshoe3_filaments = [
        BoundFilament(x1=right_gamma5, x2=left_gamma5),
        SemiInfiniteFilament(
            np.array(left_gamma5),
            va_unit,
            va_norm,
            filament_direction=1,
        ),
        SemiInfiniteFilament(
            np.array(right_gamma5),
            va_unit,
            va_norm,
            filament_direction=-1,
        ),
    ]

    horseshoes = [horseshoe1_filaments, horseshoe2_filaments, horseshoe3_filaments]

    solution = {}
    for i, horseshoe in enumerate(horseshoes):
        velocity_induced = [0, 0, 0]
        gamma = gammas[i]
        for j, filament in enumerate(horseshoe):
            logging.debug(f"j: {j},gamma: {gamma}, filament: {filament}")
            if j == 0:
                if evaluation_point[0] == filament.x1[0]:
                    tempvel = [0, 0, 0]
                else:
                    tempvel = filament.velocity_3D_bound_vortex(
                        evaluation_point, gamma, core_radius_fraction=1e-5
                    )
                logging.debug(f"BOUND tempvel: {tempvel}")
            if j == 1:
                tempvel = filament.velocity_3D_trailing_vortex_semiinfinite(
                    va_unit, evaluation_point, gamma, va_norm
                )
                logging.debug(f"SEMIINFINITE 1 tempvel: {tempvel}")
            elif j == 2:
                tempvel = filament.velocity_3D_trailing_vortex_semiinfinite(
                    va_unit, evaluation_point, gamma, va_norm
                )
                logging.debug(f"SEMIINFINITE 2 tempvel: {tempvel}")
            velocity_induced[0] += tempvel[0]
            velocity_induced[1] += tempvel[1]
            velocity_induced[2] += tempvel[2]

        solution[gammas[i]] = velocity_induced
        logging.debug(f"i: {i}")
        logging.debug(f"Analytical Solution: {analytical_solutions[gammas[i]]}")
        logging.debug(f"Numerical Solution: {solution[gammas[i]]}")

    ## creating also trailing vortices
    len_trailing = 100
    horseshoe1_filaments_trailing = [
        BoundFilament(x1=right_gamma2, x2=left_gamma2),
        BoundFilament(x1=left_gamma2, x2=left_gamma2 + len_trailing * va_unit),
        BoundFilament(
            x1=right_gamma2 + len_trailing * va_unit,
            x2=right_gamma2,
        ),
    ]
    horseshoe2_filaments_trailing = [
        BoundFilament(x1=right_gamma10, x2=left_gamma10),
        BoundFilament(x1=left_gamma10, x2=left_gamma10 + len_trailing * va_unit),
        BoundFilament(x1=right_gamma10 + len_trailing * va_unit, x2=right_gamma10),
    ]
    horseshoe3_filaments_trailing = [
        BoundFilament(x1=right_gamma5, x2=left_gamma5),
        BoundFilament(x1=left_gamma5, x2=left_gamma5 + len_trailing * va_unit),
        BoundFilament(x1=right_gamma5 + len_trailing * va_unit, x2=right_gamma5),
    ]
    horseshoes_trailing = [
        horseshoe1_filaments_trailing,
        horseshoe2_filaments_trailing,
        horseshoe3_filaments_trailing,
    ]

    logging.debug(f" ")
    logging.debug(f"TRAILING")
    solution_trailing = {}
    for i, _ in enumerate(horseshoes_trailing):
        velocity_induced = [0, 0, 0]
        gamma = gammas[i]
        for j, filament in enumerate(horseshoes_trailing[i]):
            logging.debug(f"j: {j},gamma: {gamma}, filament: {filament}")
            if j == 0:
                if evaluation_point[0] == filament.x1[0]:
                    tempvel = [0, 0, 0]
                else:
                    tempvel = filament.velocity_3D_bound_vortex(
                        evaluation_point, gamma, core_radius_fraction=1e-5
                    )
                logging.debug(f"BOUND tempvel: {tempvel}")
            if j == 1:
                tempvel = filament.velocity_3D_trailing_vortex(
                    evaluation_point, gamma, va_norm
                )
                logging.debug(f"TRAILING 1 tempvel: {tempvel}")
            elif j == 2:
                tempvel = filament.velocity_3D_trailing_vortex(
                    evaluation_point, gamma, va_norm
                )
                logging.debug(f"TRAILING 2 tempvel: {tempvel}")
            velocity_induced[0] += tempvel[0]
            velocity_induced[1] += tempvel[1]
            velocity_induced[2] += tempvel[2]

        solution_trailing[gammas[i]] = velocity_induced
        logging.debug(f"i: {i}")
        logging.debug(f"Analytical Solution: {analytical_solutions[gammas[i]]}")
        logging.debug(f"Numerical Solution: {solution_trailing[gammas[i]]}")

    return (
        solution,
        solution_trailing,
        analytical_solutions,
        gammas,
        horseshoes,
        horseshoes_trailing,
        va_unit,
    )


def test_induced_velocity_3_horsehoes():
    (
        solution,
        solution_trailing,
        analytical_solutions,
        gammas,
        horseshoe,
        horseshoes_trailing,
        va_unit,
    ) = calculate_induced_vel_3_horsheshoes()
    for i, _ in enumerate(solution):
        assert solution[gammas[i]] == pytest.approx(
            analytical_solutions[gammas[i]], abs=1e-4
        )
        assert solution_trailing[gammas[i]] == pytest.approx(
            analytical_solutions[gammas[i]], abs=1e-4
        )


def plot_line_segment(ax, segment, color, label, width: float = 3):
    ax.plot(
        [segment[0][0], segment[1][0]],
        [segment[0][1], segment[1][1]],
        [segment[0][2], segment[1][2]],
        color=color,
        label=label,
        linewidth=width,
    )
    dir = segment[1] - segment[0]
    ax.quiver(
        segment[0][0],
        segment[0][1],
        segment[0][2],
        dir[0],
        dir[1],
        dir[2],
        color=color,
    )
    return ax


def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    max_range = max([x_range, y_range, z_range])

    x_mid = np.mean(x_limits)
    y_mid = np.mean(y_limits)
    z_mid = np.mean(z_limits)

    ax.set_xlim3d([x_mid - max_range / 2, x_mid + max_range / 2])
    ax.set_ylim3d([y_mid - max_range / 2, y_mid + max_range / 2])
    ax.set_zlim3d([z_mid - max_range / 2, z_mid + max_range / 2])

    return ax


def calculate_filaments_for_plotting(filaments, chord, va_unit):
    filaments_new = []
    for i, filament in enumerate(filaments):
        logging.debug(f"filament: {filament}")
        x1 = filament.x1
        if hasattr(filament, "x2") and filament.x2 is not None:
            x2 = filament.x2
            if i == 0:  # bound
                color = "magenta"
            else:  # trailing
                color = "green"
        else:
            # For semi-infinite filaments
            x2 = x1 + 2 * chord * va_unit
            color = "orange"
            if filament.filament_direction == -1:
                x1, x2 = x2, x1
                color = "red"

        filaments_new.append([x1, x2, color])
    return filaments_new


def plot(horseshoes, va_unit):
    """
    Plots the wing panels and filaments in 3D.

    Args:
        None

    Returns:
        None
    """

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    max_chord = 1
    # Plot each panel
    for i, _ in enumerate(horseshoes):
        logging.debug(f"i: {i}, horseshoes[i]: {horseshoes[i]}")
        filaments_for_plotting = calculate_filaments_for_plotting(
            horseshoes[i], max_chord, va_unit
        )
        legends = ["Bound Vortex", "wake_1", "wake_2"]
        logging.debug(f"i: {i}, filaments_for_plotting: {filaments_for_plotting}")
        for i, _ in enumerate(filaments_for_plotting):
            x1, x2, color = filaments_for_plotting[i]
            logging.debug("Legend: %s", legends[i])
            ax = plot_line_segment(ax, [x1, x2], color, legends[i])

    # Plot the va_vector using the plot_segment
    va_vector_begin = -2 * max_chord * va_unit
    va_vector_end = va_vector_begin + 1.5 * va_unit
    ax = plot_line_segment(ax, [va_vector_begin, va_vector_end], "lightblue", "va")

    # Add legends for the first occurrence of each label
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())

    # Set equal axis limits
    ax = set_axes_equal(ax)

    # Flip the z-axis (to stick to body reference frame)
    # ax.invert_zaxis()

    # Display the plot
    plt.show()


if __name__ == "__main__":
    (
        solution,
        solution_trailing,
        analytical_solutions,
        gammas,
        horseshoes,
        horseshoes_trailing,
        va_unit,
    ) = calculate_induced_vel_3_horsheshoes()
    plot(horseshoes, va_unit)
    plot(horseshoes_trailing, va_unit)
    print(f"GAMMA order has changed due to pprint")
    print(f"Analytical Solutions")
    pprint.pprint(analytical_solutions)
    print(f"Numerical Solutions for SemiInfinite Filaments")
    pprint.pprint(solution)
    print(f"Numerical Solutions for trailing vortices")
    pprint.pprint(solution_trailing)
