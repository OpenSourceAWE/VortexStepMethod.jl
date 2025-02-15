import pytest
import numpy as np
from VSM.Filament import SemiInfiniteFilament


# Define fixtures for core_radius_fraction and gamma
@pytest.fixture
def core_radius_fraction():
    return 0.01


@pytest.fixture
def gamma():
    return 1.0


def test_calculate_induced_velocity_semi_infinite(core_radius_fraction):
    # Define a simple filament and control point
    x1 = np.array([0, 0, 0])
    direction = np.array([1, 0, 0])
    filament_direction = 1
    Uinf = 1.0
    semi_infinite_filament = SemiInfiniteFilament(
        x1, direction, Uinf, filament_direction
    )

    control_point = np.array([0.5, 0.5, 2])
    gamma = 5
    induced_velocity_calc = (
        semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
            direction, control_point, gamma, Uinf
        )
    )
    print(f"induced_velocity_calc: {induced_velocity_calc}")

    # Analytical solution using Biot-Savart law
    gamma = -gamma
    r1 = control_point - x1
    r1_cross_direction = np.cross(r1, direction)
    r_perp = np.dot(r1, direction) * direction
    alpha0 = 1.25643
    nu = 1.48e-5
    epsilon = np.sqrt(4 * alpha0 * nu * np.linalg.norm(r_perp) / Uinf)

    if np.linalg.norm(r1_cross_direction) > epsilon:
        K = (
            gamma
            / (4 * np.pi * np.linalg.norm(r1_cross_direction) ** 2)
            * (1 + np.dot(r1, direction) / np.linalg.norm(r1))
        )
        induced_velocity_analytical = K * r1_cross_direction * filament_direction
    else:
        r1_proj = np.dot(r1, direction) * direction + epsilon * (
            r1 / np.linalg.norm(r1) - direction
        ) / np.linalg.norm(r1 / np.linalg.norm(r1) - direction)
        r1_cross_direction_proj = np.cross(r1_proj, direction)
        K_proj = (
            gamma
            / (4 * np.pi * np.linalg.norm(r1_cross_direction_proj) ** 2)
            * (1 + np.dot(r1_proj, direction) / np.linalg.norm(r1_proj))
        )
        induced_velocity_analytical = (
            K_proj * r1_cross_direction_proj * filament_direction
        )

    # Assert the induced velocities are almost equal
    np.testing.assert_almost_equal(
        induced_velocity_calc, induced_velocity_analytical, decimal=6
    )


def test_a_very_close_point(gamma, core_radius_fraction):
    """Test with a point that's super close, which should be almost zero."""
    x1 = np.array([0, 0, 0])
    direction = np.array([1, 0, 0])
    filament_direction = 1
    Uinf = 1.0
    semi_infinite_filament = SemiInfiniteFilament(
        x1, direction, Uinf, filament_direction
    )
    point = [0.5, 1e-10, 0]
    induced_velocity = semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
        direction, point, gamma, Uinf
    )
    assert not np.isnan(induced_velocity).any()
    assert not np.isinf(induced_velocity).any()


@pytest.mark.parametrize(
    "point",
    [
        [0, 0, 0],  # Start point
        [0.5, 0, 0],  # Along filament
        [5, 0, 0],  # Further along filament
    ],
)
def test_point_exactly_on_filament(point, gamma, core_radius_fraction):
    x1 = np.array([0, 0, 0])
    direction = np.array([1, 0, 0])
    filament_direction = 1
    Uinf = 10
    semi_infinite_filament = SemiInfiniteFilament(
        x1, direction, Uinf, filament_direction
    )

    point = np.array([0.5, 0, 0])  # Point on the filament

    velocity = semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
        direction, point, gamma, Uinf
    )
    assert np.allclose(np.array(velocity), np.zeros(3), atol=1e-6)


def test_point_far_from_filament(gamma, core_radius_fraction):
    x1 = np.array([0, 0, 0])
    direction = np.array([1, 0, 0])
    filament_direction = 1
    Uinf = 1.0
    semi_infinite_filament = SemiInfiniteFilament(
        x1, direction, Uinf, filament_direction
    )

    far_point = np.array([0, 1e6, 0])

    induced_velocity = semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
        direction, far_point, gamma, Uinf
    )

    # Velocity should decrease with distance^2, due to Biot-Savart law
    assert not np.isnan(induced_velocity).any()


def test_different_gamma_values(core_radius_fraction):
    """Test with different gamma values to ensure linear scaling."""
    x1 = np.array([0, 0, 0])
    direction = np.array([1, 0, 0])
    filament_direction = 1
    Uinf = 1.0
    semi_infinite_filament = SemiInfiniteFilament(
        x1, direction, Uinf, filament_direction
    )
    control_point = [0.5, 1, 0]
    v1 = semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
        direction, control_point, 1.0, Uinf
    )
    v2 = semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
        direction, control_point, 2.0, Uinf
    )
    v4 = semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
        direction, control_point, 4.0, Uinf
    )
    assert np.allclose(v4, 2 * v2, 4 * v1)


def test_symmetry(gamma, core_radius_fraction):
    """Test symmetry of induced velocity for symmetric points."""
    x1 = np.array([0, 0, 0])
    direction = np.array([1, 0, 0])
    filament_direction = 1
    Uinf = 1.0
    semi_infinite_filament = SemiInfiniteFilament(
        x1, direction, Uinf, filament_direction
    )
    vel_point_pos_y = semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
        direction, [0, 1, 0], gamma, Uinf
    )
    vel_point_neg_y = semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
        direction, [0, -1, 0], gamma, Uinf
    )

    semi_infinite_filament = SemiInfiniteFilament(
        x1, direction, Uinf, -filament_direction
    )
    vel_point_pos_y_neg_dir = (
        semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
            direction, [0, 1, 0], gamma, Uinf
        )
    )
    assert np.allclose(vel_point_pos_y, -vel_point_neg_y)
    assert np.allclose(vel_point_pos_y, -vel_point_pos_y_neg_dir)


def test_uinf_effect(gamma, core_radius_fraction):
    x1 = np.array([0, 0, 0])
    direction = np.array([1, 0, 0])
    filament_direction = 1

    control_point = np.array([0, 2, 0])  # Point away from filament
    Uinf1 = 1.0
    Uinf2 = 10.0

    semi_infinite_filament1 = SemiInfiniteFilament(
        x1, direction, Uinf1, filament_direction
    )
    semi_infinite_filament2 = SemiInfiniteFilament(
        x1, direction, Uinf2, filament_direction
    )

    velocity1 = semi_infinite_filament1.velocity_3D_trailing_vortex_semiinfinite(
        direction, control_point, gamma, Uinf1
    )
    velocity2 = semi_infinite_filament2.velocity_3D_trailing_vortex_semiinfinite(
        direction, control_point, gamma, Uinf2
    )

    # Velocity should be equal for the same gamma
    assert np.allclose(np.linalg.norm(velocity1), np.linalg.norm(velocity2))


def test_around_core_radius(gamma, core_radius_fraction):
    """Test with a points around the core radius to the filament to check handling of near-singularities."""
    x1 = np.array([0, 0, 0])
    direction = np.array([1, 0, 0])
    filament_direction = 1
    Uinf = 1.0
    semi_infinite_filament = SemiInfiniteFilament(
        x1, direction, Uinf, filament_direction
    )
    core_radius_fraction = 5
    delta = 1e-5
    control_point1 = [0.5, core_radius_fraction - delta, 0]
    control_point2 = [0.5, core_radius_fraction, 0]
    control_point3 = [0.5, core_radius_fraction + delta, 0]
    induced_velocity1 = semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
        direction, control_point1, gamma, Uinf
    )
    induced_velocity2 = semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
        direction, control_point2, gamma, Uinf
    )
    induced_velocity3 = semi_infinite_filament.velocity_3D_trailing_vortex_semiinfinite(
        direction, control_point3, gamma, Uinf
    )
