"""
Internal bridge utilities for calling proprietary Siemens binaries.
"""

from __future__ import annotations

from pathlib import Path
from subprocess import Popen, PIPE
from typing import Optional, Sequence


# ---------------------------------------------------------------------------
# Locate the bin/ directory
# ---------------------------------------------------------------------------

_BIN_DIR = Path(__file__).resolve().parent.parent / "bin"

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_executable(name: str) -> Path:
    """
    Return the absolute path to an executable in the ``bin/`` folder.

    Parameters
    ----------
    name : str
        Filename of the executable or DLL, e.g. ``'SH_prototype.exe'``.

    Returns
    -------
    Path
        Absolute path to the binary.

    Raises
    ------
    FileNotFoundError
        If the file is not present in ``bin/`` — this is expected on
        machines where the proprietary Siemens tools have not been
        installed. See ``bin/README.md`` for instructions.

    Examples
    --------
    >>> from petkinetics.prototype_bridge.bridge import _get_executable
    >>> exe = _get_executable("SH_prototype.exe")
    >>> print(exe)
    /path/to/petkinetics/bin/SH_prototype.exe
    """
    exe = _BIN_DIR / "exe" / name
    if not exe.exists():
        raise FileNotFoundError(
            f"Executable '{name}' not found in {_BIN_DIR}. "
            f"See bin/README.md for instructions on obtaining proprietary tools."
        )
    return exe

def _run_executable(
    exe: Path,
    args: Sequence[str],
    check: bool = True,
) -> tuple[str, str]:
    """
    Run a binary with the given argument list and return stdout/stderr.

    Parameters
    ----------
    exe : Path
        Absolute path to the executable, typically from :func:`_get_executable`.
    args : sequence of str
        Command-line arguments to pass to the executable.
    check : bool, optional
        If ``True`` (default), raise ``RuntimeError`` when the process
        exits with a non-zero return code.

    Returns
    -------
    stdout : str
        Captured standard output from the process.
    stderr : str
        Captured standard error from the process.

    Raises
    ------
    RuntimeError
        If the process exits with a non-zero code and *check* is ``True``.
    """
    command = [str(exe)] + list(args)
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    stdout_bytes, stderr_bytes = process.communicate()

    stdout = stdout_bytes.decode("utf-8", errors="replace")
    stderr = stderr_bytes.decode("utf-8", errors="replace")

    if check and process.returncode != 0:
        raise RuntimeError(
            f"Executable '{exe.name}' exited with code {process.returncode}.\n"
            f"stderr: {stderr.strip()}"
        )

    return stdout, stderr

# ---------------------------------------------------------------------------
# Wrappers
# ---------------------------------------------------------------------------

def run_snake_VOI_extractor(
        CT_dir : str, 
        interfile_dir : str, 
        output_dir : str):
    """
    Extracts the image-derived input function of a dynamic PET scan. 

    Parameters
    ----------
    CT_dir : str
        Absolute path to CT directory containing dicom files.
    interfile_dir : str
        Absolute path to directory containing interfiles and interfile headers for each dynamic frame.
    output_dir : str
        Absolute path to where the outpu image-derived input function will be saved.

    Reference
    ----------
    Y. Tao, Z. Peng, A. Krishnan and X. S. Zhou, "Robust Learning-Based Parsing and Annotation of Medical Radiographs," 
    in IEEE Transactions on Medical Imaging, vol. 30, no. 2, pp. 338-350, Feb. 2011, doi: 10.1109/TMI.2010.2077740.
    
    """
    exe = _get_executable("getPetSnake.exe")
    stdout, stderr = _run_executable(exe, [CT_dir, interfile_dir, output_dir])

def export_raw2dcm(raw,dicom):

    """
    Converts a .raw image to a DICOM volume. 

    Parameters
    ----------
    raw : str
        Absolute path to .raw image.
    dicom : str
        Absolute path to reference dicom directory.
    """
    exe = _get_executable('rawToDICOM.exe')
    stdout, stderr = _run_executable(exe, [raw,dicom])

def simulate_tac(param_file, bif, output):

    """
    Simulates TAC with defined params file

    Parameters
    ----------
    param_file : str
        Absolute path to parameter textfile
    bif : str
        Absolute path blood input-function file
    output : str
        Absolute path output file containing simulation

    Reference
    ----------
    https://www.turkupetcentre.net/tpcclib-doc/v2/sim__3tcm_8c.html
    """

    exe = _get_executable('sim_3tcm.exe')
    vvm = 'none'
    
    _run_executable(exe, [param_file, bif, vvm, output])
