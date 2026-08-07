#!/usr/bin/env python3
import io
import sys
from contextlib import redirect_stdout, redirect_stderr


def capture_standard_output(function, *args, **kwargs):
    """Capture standard output and standard error of a given function.

    If `function` raises, the output captured up to that point is
    printed directly (to the real, un-redirected stderr) before the
    exception is re-raised, rather than being silently discarded --
    confirmed directly, from a real run, that this was NOT happening
    before: an exception from `function` propagates straight out of
    the `with redirect_stdout(...)` block, skipping past this
    function's own `return output` line entirely, so the buffered
    output (which, for a genuine SCF failure, already contains the
    full iteration trajectory up to the point of failure) was
    discarded rather than merely delayed.
    """
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()

    try:
        with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
            function(*args, **kwargs)
    except Exception:
        output = stdout_buffer.getvalue()
        errput = stderr_buffer.getvalue()
        if errput:
            output += "\n[stderr]\n" + errput
        # Printed directly here, to the REAL stderr (the with block's
        # own __exit__ has already restored sys.stderr by this point,
        # so this is not still being captured into stderr_buffer
        # itself) -- so the trajectory leading up to the failure is
        # visible even though the caller only ever sees the exception
        # itself, not this function's own (never reached) return value.
        print(output, file=sys.stderr)
        raise

    output = stdout_buffer.getvalue()
    errput = stderr_buffer.getvalue()

    if errput:
        output += "\n[stderr]\n" + errput

    return output