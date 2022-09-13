#!/usr/bin/env python
"""Examples to show xtp_binds usage."""
from pyxtp import xtp_binds
from pathlib import Path
from multiprocessing import cpu_count
import io
from contextlib import redirect_stdout


def capture_standard_output(function, *args, **kwargs):
    """Capture standard output of a given function."""
    handler = io.StringIO()
    try:
        with redirect_stdout(handler):
            function(*args, **kwargs)
    finally:
        output = handler.getvalue()
    return output
