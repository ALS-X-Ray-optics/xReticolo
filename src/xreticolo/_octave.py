"""Low-level bridge to Blazr.m running under GNU Octave.

Requires Octave on PATH (tested with 10.3.0): brew install octave

Calls are piped to a single, lazily-started, long-lived octave-cli process
rather than spawning a fresh one per call: process startup is the dominant
cost (~0.5-1s) compared to the actual RCWA computation, so reusing one
process turns a "1 second per point" sweep into "1 second, then fast".
"""

import atexit
import json
import shutil
import subprocess
import threading
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
    """Convert a Python value to an Octave/MATLAB literal for use in evaluated code."""
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(float(value))
    if isinstance(value, (list, tuple)):
        return "[" + ",".join(_to_octave_literal(v) for v in value) + "]"
    raise TypeError(f"unsupported argument type for Octave call: {type(value)!r}")


def _build_call_code(method, args, nargout):
    arg_list = ", ".join(_to_octave_literal(a) for a in args)
    return f"""
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
fflush(stdout);
"""


class OctaveSession:
    """A persistent octave-cli process that Blazr.* calls are piped to.

    Can be used as a context manager for explicit lifetime control:

        with OctaveSession() as session:
            eta = session.call("efficiency_lamellar", ...)

    Otherwise, `xreticolo.call(...)` lazily starts and reuses one shared
    session for the lifetime of the process.
    """

    def __init__(self):
        if shutil.which("octave-cli") is None:
            raise OctaveError(
                "octave-cli not found on PATH. Install Octave, e.g. `brew install octave`."
            )
        self._lock = threading.Lock()
        self._proc = subprocess.Popen(
            ["octave-cli", "--no-gui", "--norc", "--quiet"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=_RESOURCES_DIR,
        )
        self._proc.stdin.write(_SUPPRESS_WARNINGS + _SETUP_PATH)
        self._proc.stdin.flush()

    def is_alive(self):
        return self._proc.poll() is None

    def call(self, method, *args, nargout=1):
        with self._lock:
            if not self.is_alive():
                raise OctaveError("octave process is no longer running")

            self._proc.stdin.write(_build_call_code(method, args, nargout))
            self._proc.stdin.flush()

            while True:
                line = self._proc.stdout.readline()
                if line == "":
                    raise OctaveError(
                        f"octave process exited unexpectedly while running Blazr.{method}"
                    )
                if line.startswith(_ERROR_MARKER):
                    error = line[len(_ERROR_MARKER):].strip()
                    raise OctaveError(f"Blazr.{method} failed in Octave: {error}")
                if line.startswith(_RESULT_MARKER):
                    payload = line[len(_RESULT_MARKER):]
                    values = json.loads(payload)
                    outputs = [values[f"out{i}"] for i in range(1, nargout + 1)]
                    return outputs[0] if nargout == 1 else tuple(outputs)
                # anything else (warnings, stray output) is discarded

    def close(self):
        if self._proc.poll() is not None:
            return
        try:
            self._proc.stdin.close()
        except (BrokenPipeError, OSError):
            pass
        try:
            self._proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            self._proc.kill()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()


_default_session = None
_default_session_lock = threading.Lock()


def _get_default_session():
    global _default_session
    with _default_session_lock:
        if _default_session is None or not _default_session.is_alive():
            _default_session = OctaveSession()
            atexit.register(_default_session.close)
        return _default_session


def call(method, *args, nargout=1):
    """Call Blazr.<method>(*args) in Octave.

    Uses a lazily-started, process-wide persistent Octave session, so only
    the first call pays Octave's startup cost.

    Returns the single output value if nargout == 1, otherwise a tuple of
    nargout values. Numeric outputs come back as Python floats or (nested)
    lists, matching however Blazr.<method> shapes them.
    """
    return _get_default_session().call(method, *args, nargout=nargout)
