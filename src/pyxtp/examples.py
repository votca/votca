#!/usr/bin/env python
from pyxtp import xtp_binds
import os
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


def main():
    n_threads = cpu_count()
    # votca_share = Path(os.environ["VOTCASHARE"])
    # path_eanalyze = votca_share / "xtp" / "xml" / "eanalyze.xml"
    # xml_file = path_eanalyze.absolute().as_posix()
    # xtp_binds.call_calculator("eanalyze", n_threads, xml_file)

    # DFTGWBSE
    print("calling dftgwbse")
    path_dftgwbse = Path("files_examples") / "dftgwbse.xml"
    output = capture_standard_output(
        xtp_binds.call_tool, "dftgwbse", n_threads, path_dftgwbse.absolute().as_posix())


    with open("example.out", "w") as handler:
        handler.write(output)


if __name__ == "__main__":
    main()