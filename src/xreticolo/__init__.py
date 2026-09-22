"""xreticolo: Python bindings for Blazr/RETICOLO x-ray grating diffraction calculations.

Wraps Blazr.m's static methods, running under GNU Octave (no MATLAB license
required). Bundles the RETICOLO RCWA solver (Hugonin & Lalanne, CC-BY-4.0)
and a handful of tabulated material optical constants.

    import xreticolo as xr

    eta = xr.call("efficiency_lamellar", pitch_m, thickness_m, wavelength_m,
                   grazing_angle_rad, "Au")

    alpha_rad, beta_rad, theta_rad = xr.vls_trajectory(
        E_eV, E0_eV, g_lpm, p_m, q_m, cff
    )

`call()` works with any Blazr.m static method, present or future. The other
functions are thin, tested convenience wrappers around the ones used most.
"""

from ._octave import OctaveError, OctaveSession, call

__version__ = "0.1.1"


def vls_trajectory(E_eV, E0_eV, g_lpm, p_m, q_m, cff):
    """Grazing-incidence/exit angles [rad] for a VLS grating at a given energy."""
    return call("vls_trajectory", E_eV, E0_eV, g_lpm, p_m, q_m, cff, nargout=3)


def efficiency_blazedx(g_lpm, blazed_angle_deg, E_eV, grazing_angle_rad, material, order):
    """Diffraction efficiency of a blazed grating at a given order."""
    return call(
        "efficiency_blazedx",
        g_lpm, blazed_angle_deg, E_eV, grazing_angle_rad, material, order,
    )


def efficiency_lamellarx(g_lpm, thickness_m, E_eV, grazing_angle_rad, material, duty, order):
    """Diffraction efficiency of a lamellar grating at a given order."""
    return call(
        "efficiency_lamellarx",
        g_lpm, thickness_m, E_eV, grazing_angle_rad, material, duty, order,
    )


__all__ = [
    "call",
    "OctaveError",
    "OctaveSession",
    "vls_trajectory",
    "efficiency_blazedx",
    "efficiency_lamellarx",
    "__version__",
]
