# xReticolo

Python bindings for `Blazr.m`, a grating-monochromator diffraction-efficiency
tool built on top of [RETICOLO](https://zenodo.org/records/14631951) (RCWA
solver by J.P. Hugonin and P. Lalanne). Calculations run in
[GNU Octave](https://octave.org/) — no MATLAB license required.

## Requirements

- Python >= 3.9
- [Octave](https://octave.org/) on `PATH` (tested with 10.3.0):
  `brew install octave` (macOS) / `apt install octave` (Debian/Ubuntu)

## Install

```bash
pip install xreticolo
```

## Quickstart

```python
import xreticolo as xr

# blazed grating, +1 order efficiency at 806 eV
alpha_rad, beta_rad, theta_rad = xr.vls_trajectory(
    E_eV=806, E0_eV=806, g_lpm=287.44e3, p_m=24.301, q_m=7.573, cff=1.577
)
eta = xr.efficiency_blazedx(
    g_lpm=287.44e3, blazed_angle_deg=0.40, E_eV=806,
    grazing_angle_rad=1.5707963267948966 - alpha_rad, material="Au", order=1,
)
print(eta)
```

`xr.call(method, *args, nargout=1)` invokes any `Blazr.<method>` static
method directly, so it keeps working as `Blazr.m` gains new methods:

```python
eta = xr.call("efficiency_lamellar", pitch_m, thickness_m, wavelength_m,
              grazing_angle_rad, "Au")
```

See [`examples/`](examples/) for full worked notebooks (blazed and lamellar
gratings, efficiency vs. photon energy).

## What's bundled

- `Blazr.m` — the grating efficiency/trajectory wrapper.
- The RETICOLO V7 RCWA engine (`reticolo_allege`), by J.P. Hugonin and
  P. Lalanne, CC-BY-4.0 — see
  [arXiv:2101.00901](https://arxiv.org/abs/2101.00901) and
  [zenodo.org/records/14631951](https://zenodo.org/records/14631951).
  If you use this package for published work, please cite that paper.
- A handful of tabulated optical constants (Au, Ni, Rh, SiC, C) used by
  `Blazr.m`'s efficiency functions.

## Development

```bash
pip install -e ".[examples,test]"
pytest
```

## License

BSD-3-Clause (see `LICENSE`) for the Python bindings and `Blazr.m`. The
bundled RETICOLO engine is distributed under its original CC-BY-4.0 license
from the authors.
