"""Low-level bridge to Blazr.m running under GNU Octave.

Requires Octave on PATH (tested with 10.3.0): brew install octave
"""

import json
import shutil
import subprocess
from pathlib import Path

_PACKAGE_DIR = Path(__file__).resolve().parent
_RESOURCES_DIR = _PACKAGE_DIR / "_resources"
_ENGINE_DIR = _RESOURCES_DIR / "RETICOLO V7"
_MATERIAL_DIR = _RESOURCES_DIR / "material_data"

# RETICOLO is old MATLAB code (French comments, & instead of &&, etc.);
# these warnings are cosmetic noise, not real problems.
_SUPPRESS_WARNINGS = """
warning('off','Octave:possible-matlab-short-circuit-operator');
warning('off','Octave:legacy-function');
warning('off','Octave:num-to-str');
warning('off','octave:get_input:invalid_utf8');
warning('off','Octave:data-file-in-path');
"""

_SETUP_PATH = f"""
addpath('{_RESOURCES_DIR}');
addpath('{_MATERIAL_DIR}');
addpath(genpath('{_ENGINE_DIR}'));
"""

_RESULT_MARKER = "XRETICOLO_JSON:"
_ERROR_MARKER = "XRETICOLO_ERROR:"


class OctaveError(RuntimeError):
    """Raised when a Blazr.m call fails inside Octave."""


def _to_octave_literal(value):
    """Convert a Python value to an Octave/MATLAB literal for use in --eval'd code."""
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(float(value))
    if isinstance(value, (list, tuple)):
        return "[" + ",".join(_to_octave_literal(v) for v in value) + "]"
    raise TypeError(f"unsupported argument type for Octave call: {type(value)!r}")


def _extract(stdout, marker):
    for line in stdout.splitlines():
        if line.startswith(marker):
            return line[len(marker):]
    return None


def call(method, *args, nargout=1):
    """Call Blazr.<method>(*args) in Octave.

    Returns the single output value if nargout == 1, otherwise a tuple of
    nargout values. Numeric outputs come back as Python floats or (nested)
    lists, matching however Blazr.<method> shapes them.
    """
    if shutil.which("octave-cli") is None:
        raise OctaveError(
            "octave-cli not found on PATH. Install Octave, e.g. `brew install octave`."
        )

    arg_list = ", ".join(_to_octave_literal(a) for a in args)
    code = f"""
try
  call_args = {{{arg_list}}};
  outs = cell(1, {nargout});
  [outs{{:}}] = Blazr.{method}(call_args{{:}});
  result = struct();
  for i_out = 1:{nargout}
    result.(sprintf('out%d', i_out)) = outs{{i_out}};
  end
  printf('{_RESULT_MARKER}%s\\n', jsonencode(result));
catch err
  printf('{_ERROR_MARKER}%s\\n', strrep(err.message, "\\n", " "));
end
"""
    result = subprocess.run(
        ["octave-cli", "--no-gui", "--eval", _SUPPRESS_WARNINGS + _SETUP_PATH + code],
        cwd=_RESOURCES_DIR,
        capture_output=True,
        text=True,
    )

    error = _extract(result.stdout, _ERROR_MARKER)
    if error is not None:
        raise OctaveError(f"Blazr.{method} failed in Octave: {error.strip()}")

    payload = _extract(result.stdout, _RESULT_MARKER)
    if payload is None:
        raise OctaveError(
            f"no result from Blazr.{method} (octave exit {result.returncode})\n"
            f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
        )

    values = json.loads(payload)
    outputs = [values[f"out{i}"] for i in range(1, nargout + 1)]
    return outputs[0] if nargout == 1 else tuple(outputs)
