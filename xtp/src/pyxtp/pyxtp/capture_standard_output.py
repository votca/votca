#!/usr/bin/env python3
import io
from contextlib import redirect_stdout, redirect_stderr


def capture_standard_output(function, *args, **kwargs):
    """Capture standard output and standard error of a given function."""
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()

    with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
        function(*args, **kwargs)

    output = stdout_buffer.getvalue()
    errput = stderr_buffer.getvalue()

    if errput:
        output += "\n[stderr]\n" + errput

    return output