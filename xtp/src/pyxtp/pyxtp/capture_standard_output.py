#!/usr/bin/env python
"""Examples to show xtp_binds usage."""
import io
from contextlib import redirect_stdout


def capture_standard_output(function, *args, **kwargs):
    """Capture standard output of a given function."""
    with redirect_stdout(io.StringIO()) as stdout:
        function(*args, **kwargs)

    return stdout.getvalue()
