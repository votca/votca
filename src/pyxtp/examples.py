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


def run_mapchecker(nThreads: int, path_examples: Path) -> None:
    """Run the mapchecker calculator."""
    # mapchecker
    print("Calling mapchecker")
    # Look for the xml file to use as input
    xml_file = (path_examples / "mapchecker.xml").absolute().as_posix()
    state_file = (path_examples / "state.hdf5").absolute().as_posix()
    # Input to call the calculator
    dict_inp = {
        "xmlFile": xml_file,
        "stateFile": state_file,
        "nThreads": str(nThreads),
        "nFrames": '1',
        "firstFrame": '0',
        "save": '1'  # xtp_binds casts it to bool False
    }
    xtp_binds.call_calculator("mapchecker", dict_inp)

    # Check that the files exists
    expected = {"md_segments_step_0.pdb", "mp_segments_e_step_0.pdb",
                "mp_segments_h_step_0.pdb", "qm_segments_n_step_0.pdb"}
    pdb_files = set(p.as_posix() for p in Path(".").glob("*_step_0.pdb"))
    assert all(x in expected for x in pdb_files)


def run_dftgwbse(nThreads: int, path_examples: Path) -> None:
    """Run the DFTGWBSE qmtool."""
    print("calling dftgwbse")
    # Look for the xml input
    path_dftgwbse = (path_examples / "dftgwbse.xml").absolute().as_posix()
    # Call the tool and capture the standard output
    output = capture_standard_output(
        xtp_binds.call_tool, "dftgwbse", nThreads, path_dftgwbse)

    with open("example.out", "w") as handler:
        handler.write(output)


def run_examples():
    """Call the xtp_binds interface."""
    nThreads = cpu_count()
    path_examples = Path("files_examples")

    run_mapchecker(nThreads, path_examples)
    run_dftgwbse(nThreads, path_examples)


if __name__ == "__main__":
    run_examples()
