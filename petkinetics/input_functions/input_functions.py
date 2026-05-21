from __future__ import annotations

from pathlib import Path
import os
import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from petkinetics.prototype_bridge import bridge


# ---------------------------------------------------------------------------
# Locate the data/ directory
# ---------------------------------------------------------------------------

_DATA_DIR = Path(__file__).resolve().parent / "data"

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _write_tac(time, activity, dst):
    with open (dst,'a') as proc_seqf:                   # Create and open txt file with the given path
        for t, a in zip(time, activity):       # For every value in the given x- and y-values:
            proc_seqf.write("{}\t{}\n".format(t, a[0]))    # Write the x- and y-values seperated by a tab followed by a new line

def _plot_sPBIF(time, activity, dst):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5),gridspec_kw={'width_ratios': [6,4]})

    ax1.plot(time, np.divide(activity, 1000), color='k')
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Activity [kBq]')

    ax2.plot(time, np.divide(activity, 1000), color='k')
    ax2.set_xlabel('Time [s]')
    ax2.set_xlim(0, 300)  # Zoom in on the first 100 seconds

    fig.suptitle('sPBIF')
    plt.tight_layout()
    plt.savefig(rf'{dst}\sPBIF.png', dpi = 96.0)

# ---------------------------------------------------------------------------
# Public API
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
    bridge.run_snake_VOI_extractor(
        CT_dir, 
        interfile_dir, 
        output_dir)

def gen_sPBIF(
        time: NDArray[np.float64],
        activity: NDArray[np.float64],
        dst: str = None,
        plot: bool = False,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Generate a scaled population-based input function (sPBIF).

    Loads the population-based input function from Dias et al. (2022),
    then scales it to match the measured tissue or blood activity by
    computing a single scaling factor from the ratio of mean activities
    over the supplied time points.

    Parameters
    ----------
    time : NDArray[np.float64]
        Mid-frame times in minutes, shape ``(n_frames,)``. Used to look
        up matching time points in the population-based input function
        and to compute the mean activity for scaling. Must be a subset
        of the time points defined in the Dias et al. data file.
    activity : NDArray[np.float64]
        Measured activity concentrations (e.g. kBq/mL) at each time
        point in *time*, shape ``(n_frames,)``. Typically an image-derived
        or blood sample TAC used as the scaling reference.
    dst : str, optional
        File path for saving the sPBIF to disk. If ``None`` (default),
        no file is written. See :func:`_save_sPBIF` for the output format.
    plot : bool, optional
        If ``True``, generate a plot of the scaled PBIF via
        :func:`_plot_sPBIF`. Requires *dst* to be set if the plot is to
        be saved. Exceptions during plotting are caught and printed
        without interrupting execution. Default is ``False``.

    Returns
    -------
    PB_Time : NDArray[np.float64]
        Full time axis of the population-based input function in minutes,
        shape ``(n_pbif_frames,)``. Covers the complete duration defined
        in the Dias et al. data file, not just the frames in *time*.
    sPBIF_activity : NDArray[np.float64]
        Scaled activity values corresponding to *PB_Time*, length
        ``n_pbif_frames``.

    Notes
    -----
    The scaling factor is defined as:

    .. math::

        s = \\frac{\\bar{C}_{\\text{measured}}}{\\bar{C}_{\\text{PBIF}}}

    where both means are computed only over the time points supplied in
    *time*. The factor is then applied uniformly across the full PBIF
    time axis.

    The population-based input function is loaded from
    ``input_functions/data/dias2022_PBIF.txt``, which is included with
    the package. Do not modify this file.

    References
    ----------
    Dias AH, Smith AM, Shah V, Pigg D, Gormsen LC, Munk OL. Clinical 
    validation of a population-based input function for 20-min dynamic 
    whole-body 18F-FDG multiparametric PET imaging. EJNMMI Phys. 2022 
    Sep 8;9(1):60. doi: 10.1186/s40658-022-00490-y.
    """
    # ---- Load PBIF ------------------------------------------------
    PB_Time, PB_activity = [], []
    with open(_DATA_DIR / "dias2022_PBIF.txt", "r") as f:
        lines=f.readlines(); lines = lines[1:len(lines)]
        for x in lines:
            PB_Time.append(float(x.split('\t')[0]))
            PB_activity.append(float(x.split('\t')[1]))

    PB_Time = np.array(PB_Time); 
    PB_activity = np.array(PB_activity); 
    
    # ---- Scale PBIF ----------------------------------------------

    activity_sPBIF=[]
    [activity_sPBIF.append(PB_activity[PB_Time==i]) for i in time]
    
    # Calculate averages
    avg_act = sum(activity)/len(time) 
    avg_act_PBIF = sum(activity_sPBIF)/len(activity_sPBIF)

    # Calculate scaling factor for the scaled PBIF (sPBIF) 
    scaling_factor = avg_act/avg_act_PBIF  

    # Multiply each element with the scaling factor                            
    sPBIF_activity = np.array([scaling_factor * element for element in PB_activity])
    
    if dst:
        _write_tac(PB_Time, sPBIF_activity, dst)
    if plot:
        try:
            _plot_sPBIF(PB_Time, sPBIF_activity, dst)
        except Exception as e:
            print({e})

    return PB_Time, sPBIF_activity
