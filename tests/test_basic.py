import math
import shutil

import pytest

xr = pytest.importorskip("xreticolo")

pytestmark = pytest.mark.skipif(
    shutil.which("octave-cli") is None, reason="octave-cli not on PATH"
)


def test_vls_trajectory_and_blazed_efficiency():
    alpha_rad, beta_rad, theta_rad = xr.vls_trajectory(
        E_eV=806, E0_eV=806, g_lpm=287.44e3, p_m=24.301, q_m=7.573, cff=1.577
    )
    assert 1.0 < alpha_rad < math.pi / 2

    eta = xr.efficiency_blazedx(
        g_lpm=287.44e3, blazed_angle_deg=0.40, E_eV=806,
        grazing_angle_rad=math.pi / 2 - alpha_rad, material="Au", order=1,
    )
    assert 0 < eta < 1


def test_efficiency_lamellarx():
    alpha_rad, _, _ = xr.vls_trajectory(
        E_eV=806, E0_eV=806, g_lpm=179e3, p_m=24.301, q_m=7.573, cff=1.632
    )
    eta = xr.efficiency_lamellarx(
        g_lpm=179e3, thickness_m=19e-9, E_eV=806,
        grazing_angle_rad=math.pi / 2 - alpha_rad, material="Au", duty=0.38, order=1,
    )
    assert 0 < eta < 1


def test_call_generic_matches_convenience_wrapper():
    direct = xr.call(
        "efficiency_blazedx", 287.44e3, 0.40, 806, 1.4, "Au", 1,
    )
    convenience = xr.efficiency_blazedx(287.44e3, 0.40, 806, 1.4, "Au", 1)
    assert direct == convenience
