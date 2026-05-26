import os
from subprocess import Popen
import pandas as pd
from pathlib import Path
import uuid
import numpy as np
import sys
import matplotlib.pyplot as plt
import subprocess
import math as ma
from contextlib import contextmanager
from itertools import islice
from petkinetics.model_fitting import model_fitting as mf 
from petkinetics.prototype_bridge import bridge

leaveoutDP1 = True
saveHighBiasFits = True

SIM_3TCM = r"D:\PhD_project_JRM\Software\TPC-software\bin\sim_3tcm.exe"

@contextmanager
def _suppress_c_stdout():
    devnull = os.open(os.devnull, os.O_WRONLY)
    old_stdout = os.dup(1)
    # old_stderr = os.dup(2)   # ← leave stderr open for now
    try:
        os.dup2(devnull, 1)
        # os.dup2(devnull, 2)  # ← leave stderr open for now
        yield
    finally:
        os.dup2(old_stdout, 1)
        os.close(devnull)
        os.close(old_stdout)

# ─────────────────────────────────────────────────────────────────────────────
# TRUE CALCULATION FUNCTIONS
#
# These compute the ground-truth distribution volume (Vt) from kinetic
# parameters. Two variants exist depending on which parameter pair is varied.
# ─────────────────────────────────────────────────────────────────────────────

def _logan_vt(sim_params: dict, vb: float) -> float:
    """
    Compute the theoretical total volume of distribution (Vt) for the Logan plot.

    Derived analytically from a two-tissue compartment model:

        Vt = (1 - vb/100) * (K1/k2) * (1 + k3/k4) + vb/100

    Used to generate ground-truth Vt values for simulated TACs, against
    which Logan-estimated Vt can be validated.

    Parameters
    ----------
    sim_params : dict
        Kinetic rate constants with keys 'K1', 'k2', 'k3', 'k4',
        in units of min⁻¹ (K1 in ml/cm³/min).
    vb : float
        Blood volume fraction as a percentage (0–100).

    Returns
    -------
    float
        Theoretical Vt in ml/cm³, rounded to 4 decimal places.
    """
    K1 = sim_params['K1']
    k2 = sim_params['k2']
    k3 = sim_params['k3']
    k4 = sim_params['k4']
    return round((1 - vb / 100) * (K1 / k2) * (1 + k3 / k4) + vb / 100, 4)


def _patlak_ki(sim_params: dict, vb: float) -> float:
    """
    Compute the theoretical influx rate constant (Ki) for the Patlak plot.

    Derived analytically from a two-tissue compartment model with irreversible
    trapping (k4 = 0):

        Ki = (1 - vb/100) * (K1*k3) / (k2 + k3) + vb/100

    Used to generate ground-truth Ki values for simulated TACs, against
    which Patlak-estimated Ki can be validated.

    Parameters
    ----------
    sim_params : dict
        Kinetic rate constants with keys 'K1', 'k2', 'k3',
        in units of min⁻¹ (K1 in ml/cm³/min). k4 is assumed zero
        and should not be included.
    vb : float
        Blood volume fraction as a percentage (0–100).

    Returns
    -------
    float
        Theoretical Ki in ml/cm³/min, rounded to 4 decimal places.
    """
    K1 = sim_params['K1']
    k2 = sim_params['k2']
    k3 = sim_params['k3']
    return round((1 - vb / 100) * ((K1 * k3) / (k2 + k3)) + vb / 100, 4)

# Dispatch tables for the true value calculations. Maps the user-facing method name to the correct function.
METHOD_DISPATCHERS = {
        'calcTrueVt':    _logan_vt,
        'calcTrueKi':    _patlak_ki,
    }

# ─────────────────────────────────────────────────────────────────────────────
# MODEL FITTING FUNCTIONS
#
# These compute the fitted  kinetic parameters.
# ─────────────────────────────────────────────────────────────────────────────

def fit_model(model_type, tac_path, bif_path, get_tac = False):
    slope, ic, error, model, frameTime, tac = MODEL_DISPATCHERS[model_type](
                tac_path, bif_path, figFlag = False)
    if get_tac:
        return slope, ic, error, model, frameTime, tac
    else:
        return slope, ic, error, model

def dLogan_tStar(tac_path, bif_path, figFlag = False):
    return mf.dLoganRB_tStar(tac_path, bif_path, k1k2Flag = True, figFlag=figFlag)

def dLogan_tZero(tac_path, bif_path, figFlag = False):
    return mf.dLoganRB_tZero(tac_path, bif_path, k1k2Flag = False, figFlag=figFlag)

def dLogan_tStar_alphaBeta(tac_path, bif_path, figFlag = False):
    return mf.dLoganRB_tStar(tac_path, bif_path,  k1k2Flag = False, figFlag=figFlag)

def dLogan_tZero_K1k2(tac_path, bif_path, figFlag = False):
    return mf.dLoganRB_tZero(tac_path, bif_path, k1k2Flag = True, figFlag=figFlag)

def Patlak_tStar_Linear(tac_path, bif_path, figFlag = False):
    return mf.Patlak_tStar(tac_path, bif_path, figFlag=figFlag, linear=True)

def Patlak_tStar_NonLinear(tac_path, bif_path, figFlag = False):
    return mf.Patlak_tStar(tac_path, bif_path, figFlag=figFlag, linear=False)

def Patlak_tZero(tac_path, bif_path, figFlag = False):
    return mf.Patlak_tZero(tac_path, bif_path, figFlag=figFlag)

# The dispatch table maps the user-facing model_type string to the correct fitting function.
MODEL_DISPATCHERS = {
    "dLogan_tStar":             dLogan_tStar,
    "dLogan_tZero":             dLogan_tZero,
    "dLogan_tStar_alphaBeta":   dLogan_tStar_alphaBeta,
    "dLogan_tZero_K1k2":        dLogan_tZero_K1k2,
    "Patlak_tStar_linear":      Patlak_tStar_Linear,
    "Patlak_tStar_nonlinear":   Patlak_tStar_NonLinear,
    "Patlak_tZero":             Patlak_tZero,
}


# ─────────────────────────────────────────────────────────────────────────────
# TAC SIMULATION
# ─────────────────────────────────────────────────────────────────────────────


def _convert_sec2min(input_file, output_file):
    # Read the file
    with open(input_file, "r") as f:
        lines = f.readlines()

    converted_lines = []

    for line in lines:
        # Skip empty lines
        if not line.strip():
            continue

        # Split by arbitrary whitespace (tab or spaces)
        parts = line.split()

        # Expect two columns: time and activity
        if len(parts) < 2:
            continue

        # Read values
        time_sec = float(parts[0])
        activity = parts[1]

        # Convert seconds to minutes
        time_min = time_sec/60

        # Write with tab separation
        converted_lines.append(f"{time_min:.6f}\t{activity}\n")

        # Save new file
    with open(output_file, "w") as f:
        f.writelines(converted_lines)

def _convert_min2sec(input_file, output_file):
    # Read the file
    with open(input_file, "r") as f:
        lines = f.readlines()

    converted_lines = []

    for line in lines:
        # Skip empty lines
        if not line.strip():
            continue

        # Split by arbitrary whitespace (tab or spaces)
        parts = line.split()

        # Expect two columns: time and activity
        if len(parts) < 2:
            continue

        # Read values
        time_sec = float(parts[0])
        activity = parts[1]

        # Convert seconds to minutes
        time_min = round(time_sec*60)

        # Write with tab separation
        converted_lines.append(f"{time_min:.6f}\t{activity}\n")
        
        # Save new file
    with open(output_file, "w") as f:
        f.writelines(converted_lines)

def simulateTac(bif: str, simfile: str, params: dict, temp_dir: str, vb) -> None:
    """
    Simulate a TAC using the external 3TCM binary.

    Parameters
    ----------
    bif : str
        Path to the blood input function file.
    simfile : str
        Path where the output TAC will be written.
    params : dict
        Compartment model parameters. Valid keys: K1, k2, k3, k4, k5, k6, Vb.
        Any key not provided defaults to the values below.
    """
    # Start from default values, then overwrite with whatever was provided.
    defaults = dict(K1=0.8, k2=0.98, k3=0.0, k4=0.0, k5=0.0, k6=0.0, Vb=vb)
    defaults.update(params)
    p = defaults

    # uuid4() generates a random unique identifier — this ensures each
    # parallel worker writes to a different param file and avoids overwrites.
    param_file = os.path.join(temp_dir, f"params_{uuid.uuid4().hex}.txt")

    header = (
        "# model := SER3TCM\n"
        "Parameters\tK1[mL/(mL*min)]\tk2[1/min]\tk3[1/min]\t"
        "k4[1/min]\tk5[1/min]\tk6[1/min]\tVb[%]\n"
    )
    row = (
        f"Tac\t{p['K1']}\t{p['k2']}\t{p['k3']}\t"
        f"{p['k4']}\t{p['k5']}\t{p['k6']}\t{p['Vb']}"
    )

    with open(param_file, 'w') as f:
        f.write(header + row)

    # Remove any existing output file to avoid the binary appending to it.
    if os.path.exists(simfile):
        os.remove(simfile)

    # Run the external binary. Popen().communicate() blocks until it finishes.
    bridge.simulate_tac(param_file, bif, simfile)

    # The binary writes a redundant header line — strip it out.
    with open(simfile, "r") as f:
        lines = [line for line in f if line.strip() != "time[unknown]\tTac[unknown]"]
    # lines = [line.replace("\t"," ") for line in lines if "\t" in line]
    with open(simfile, "w") as f:
        f.writelines(lines)

    _convert_min2sec(simfile,simfile)

    # Clean up the temporary parameter file.
    os.remove(param_file)

# ─────────────────────────────────────────────────────────────────────────────
# TAC FILE LOADER AND SAVER
# ─────────────────────────────────────────────────────────────────────────────

def _load_tac(simfile: str):
    """
    Load a two-column tab-separated TAC file.

    Returns
    -------
    time : np.ndarray
    activity : np.ndarray
    """
    data = pd.read_csv(simfile, sep='\t', header=None, names=['time', 'activity'])
    return data['time'].values, data['activity'].values

def _write_tac(path: str, time: np.ndarray, activity: np.ndarray) -> None:
    """Write a two-column TAC file (tab-separated, no header)."""
    with open(path, 'w') as f:
        for t, a in zip(time, activity):
            f.write(f"{t}\t{a}\n")

# ─────────────────────────────────────────────────────────────────────────────
# IN-MEMORY NOISE HELPER
# ─────────────────────────────────────────────────────────────────────────────

def _apply_relative_noise(
    activity: np.ndarray,
    fraction: float,
    seed: int = 42,
) -> np.ndarray:
    """
    Add Gaussian noise whose standard deviation is proportional to the local
    signal amplitude (relative/SNR-preserving noise model).

    TAC_noisy(t) = TAC(t) + N(0, fraction * |TAC(t)|)

    Parameters
    ----------
    activity : Clean TAC values.
    fraction : Noise level as a fraction of signal amplitude (e.g. 0.05 = 5%).
    seed     : Random seed for reproducibility.
    """
    rng       = np.random.default_rng(seed)
    # Use absolute value of signal as the local scale to avoid negative std.
    local_std = fraction * np.abs(activity)
    noise     = rng.normal(loc=0.0, scale=local_std)
    return activity + noise


# ─────────────────────────────────────────────────────────────────────────────
# WORKER FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

def _worker(task: dict) -> dict | None:
    """
    Simulate once tac; fit across all (noise_level, tstart) combinations.

    """

    uid     = uuid.uuid4().hex
    simfile = str(task['simfile_template']).replace('.txt', f"_{uid}.txt")

    # ── 1. Simulate once ──────────────────────────────────────────────────────
    try:
        simulateTac(task['framebased_bif'], simfile, task['sim_params'], task['temp_dir'], task['vb'])
    except Exception as e:
        import traceback
        print(f"[Worker] Simulation failed at i={task['i']}, j={task['j']}, "
            f"l={task['l']}, m={task['m']}: {e}", file=sys.__stderr__)
        traceback.print_exc(file=sys.__stderr__)
        return None

    # ── 2. Load noiseless TAC ─────────────────────────────────────────────────
    time, activity_clean = _load_tac(simfile)

    # ── 3. Ground-truth values ────────────────────────────────────────────────────

    true_values = dict()
    if task['methods']['Logan']:
        Vt_true = METHOD_DISPATCHERS['calcTrueVt'](task['sim_params'],task['vb'])
        true_values.update(Logan_Vt = Vt_true)
    if task['methods']['Patlak']:
        Ki_true = METHOD_DISPATCHERS['calcTrueKi'](task['sim_params'],task['vb'])
        true_values.update(Patlak_Ki = Ki_true)

    n_dp  = task['datapoints']

    noise_tstart_results = {}
    noise_tstart_tac = {}

    # ── 4. Noise loop ─────────────────────────────────────────────────────────
    for noise_idx, noise_level in enumerate(task['noise_list']):

        if noise_level > 0:
            activity = _apply_relative_noise(
                activity_clean, noise_level, seed = int(42 + noise_idx * 10000 + task['i'] * 100 + task['j'])
            )
        else:
            activity = activity_clean

        # Unique noisy TAC file per worker × noise level — no collision possible
        noisy_file = simfile.replace('.txt', f"_ni{noise_idx}.txt")
        _write_tac(noisy_file, time, activity)

        tstart_results = {}
        tstart_tac = {}

        # ── 5. tstart loop ────────────────────────────────────────────────────
        for tstart_sec in task['tstart_seconds']:

            # Build a unique tacL8 path for this worker × tstart combination
            tac_base  = task['tacL8_base']
            tac_DP  = f"{tac_base}_ts{tstart_sec}_{uid}_ni{noise_idx}.txt"

            start_idx = int(np.searchsorted(time, tstart_sec, side='left'))
            start_idx = min(start_idx, len(time) - n_dp)

            _write_tac(tac_DP, 
                       time[start_idx : start_idx + n_dp], 
                       activity[start_idx : start_idx + n_dp])

            tac_saved = 0

            model_collection = []
            if task['methods']['Logan']: model_collection += task['methods']['Logan']
            if task['methods']['Patlak']: model_collection += task['methods']['Patlak']

            model_collection_result = {}
            for model_type in model_collection:

                with _suppress_c_stdout():
                    if tac_saved == 0:
                        slope, ic, error, model, frameTime, tac = fit_model(
                            model_type, tac_DP, task['interpolated_bif'], get_tac=True)
                        tac_saved = 1
                    else:
                        slope, ic, error, model = fit_model(
                            model_type, tac_DP, task['interpolated_bif'], get_tac=False)
                        
                model_collection_result[model_type] = dict(slope=slope, ic=ic, error=error, model=model)
            
            tstart_tac[tstart_sec]     = dict(time = frameTime, tac = tac)
            tstart_results[tstart_sec]  = model_collection_result

            if os.path.exists(tac_DP):
                os.remove(tac_DP) 
        
        if os.path.exists(noisy_file):
            os.remove(noisy_file)   

        noise_tstart_results[noise_level] = tstart_results
        noise_tstart_tac[noise_level] = tstart_tac
    os.remove(simfile)

    return dict(
        i                    = task['i'],
        j                    = task['j'],
        l                    = task['l'],
        m                    = task['m'],
        true_values          = true_values,
        noise_tstart_results = noise_tstart_results,
        noise_tstart_tac     = noise_tstart_tac,
    )


# ─────────────────────────────────────────────────────────────────────────────
#  MAIN ORCHESTRATOR
# ─────────────────────────────────────────────────────────────────────────────

def _memmap_path(orgDir, vb, nl, ts, kind):
    """Canonical path for one memmap .npy file."""
    return Path(orgDir) / (
        f"sim_{kind}_vb{vb}_noiseLvl{nl}_tstart{int(ts)}.npy"
    )

def _tac_path(orgDir, vb, nl, ts):
    """Canonical path for one memmap .npy file."""
    return Path(orgDir) / (
        f"sim_tac_vb{vb}_noiseLvl{nl}_tstart{int(ts)}.npy"
    )

def _sentinel_path(orgDir, vb):
    """Path for the boolean completion sentinel array."""
    return Path(orgDir) / f"_sentinel_vb{vb}.npy"

def _result_path(orgDir, vb):
    return Path(orgDir) / f"sim_result_vb{vb}.npy"

def _params_path(orgDir):
    return Path(orgDir) / f"sim_params.npy"

def _sim_path(temp, K1, k2, k3, k4):
    return Path(temp) / f"sim_K1_{K1}_k2_{k2}_k3_{k3}_k4_{k4}.txt"

def _SS_path(temp, K1, k2, k3, k4):
    return Path(temp) / f"sim_Steady-State_K1_{K1}_k2_{k2}_k3_{k3}_k4_{k4}.txt"

def simulate_and_fit(
    K1_list,k2_list,k3_list,k4_list,
    vb,
    methods,
    framebased_bif, interpolated_bif,
    tstart_seconds,
    orgDir,
    temp_dir,
    datapoints = 8,
    n_workers  = 9,
    noise_list = [0.0]):


    N1, N2, N3, N4  = len(K1_list), len(k2_list), len(k3_list), len(k4_list)
    orgDir          = Path(orgDir)

    # ── Open or create memory-mapped arrays ───────────────────────────────────

    fit_p  = _result_path(orgDir, vb)
    sen_p = _sentinel_path(orgDir, vb)
    param_p = _params_path(orgDir)

    mode_vt  = 'r+' if fit_p.exists()  else 'w+'
    mode_param = 'r+' if param_p.exists() else 'w+'
    mode_sen = 'r+' if sen_p.exists() else 'w+'

    mm_true  = np.lib.format.open_memmap(str(fit_p),  mode=mode_vt,
                                        dtype=float, shape=(N1, N2, N3, N4, 2))
    mm_params  = np.lib.format.open_memmap(str(param_p),  mode=mode_param,
                                        dtype=float, shape=(N1, N2, N3, N4, 4))
    mm_sen = np.lib.format.open_memmap(str(sen_p), mode=mode_sen,
                                        dtype=bool,  shape=(N1, N2, N3, N4))

    if mm_true.mode == 'w+':
        mm_true[:] = np.nan

    # One memmap per (noise_level, tstart, kind) combination
    model_collection = []
    if methods['Logan']: model_collection += methods['Logan']
    if methods['Patlak']: model_collection += methods['Patlak']

    mm = {}
    mm_tac = {}
    for nl in noise_list:
        for ts in tstart_seconds:

            tac_p = _tac_path(orgDir, vb, nl, ts)
            mode = 'r+' if tac_p.exists() else 'w+'
            arr_tac = np.lib.format.open_memmap(
                str(tac_p), mode=mode, dtype=float,
                shape=(N1, N2, N3, N4, 2, datapoints)
            )
            if mode == 'w+':
                arr_tac[:] = np.nan
            mm_tac[(nl, ts)] = arr_tac

            for model_type in model_collection:
                p    = _memmap_path(orgDir, vb, nl, ts, model_type)
                mode = 'r+' if p.exists() else 'w+'
                arr  = np.lib.format.open_memmap(
                    str(p), mode=mode, dtype=float,
                    shape=(N1, N2, N3, N4, 3 + datapoints)) # 4 for both fitted slope, ic, error and model values

                if mode == 'w+':
                    arr[:] = np.nan

                mm[(nl, ts, model_type)] = arr


    # ── Build task list, skipping already-completed cells ─────────────────────
    
    tasks = []
    skipped = 0
    for i, K1 in enumerate(K1_list):

        for j, k2 in enumerate(k2_list):
            for l, k3 in enumerate(k3_list):
                for m, k4 in enumerate(k4_list):
                    if mm_sen[i, j, l, m]:       # already saved on a previous run
                        skipped += 1
                        continue
                    tasks.append(dict(
                    i=i, j=j, l=l, m=m,
                    vb=vb,
                    sim_params=dict(K1=K1, k2=k2, k3=k3, k4=k4),
                    methods=methods,
                    framebased_bif=framebased_bif,
                    interpolated_bif=interpolated_bif,
                    simfile_template=_sim_path(temp_dir, K1, k2, k3, k4),
                    noise_list=noise_list,
                    tstart_seconds=tstart_seconds,
                    tacL8_base=_SS_path(temp_dir, K1, k2, k3, k4),
                    datapoints=datapoints,
                    temp_dir=temp_dir
                ))



    if skipped:
        print(f"    [RESUME] {skipped}/{N1*N2*N3*N4} grid points already done, "
              f"submitting {len(tasks)} remaining.")

    if not tasks:
        print(f"    [SKIP] All {N1*N2*N3*N4} grid points complete.")
        return


    # ── Parallel execution with immediate incremental writes ──────────────────

    from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED

    completed  = 0
    total      = len(tasks)
    in_flight_tasks = {}   # future → task
    MAX_QUEUED = n_workers * 4

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        task_iter = iter(tasks)

        # Seed the initial window
        print(f"    [SUBMIT] Seeding initial {min(MAX_QUEUED, total)} tasks...")
        for t in islice(task_iter, MAX_QUEUED):
            fut = executor.submit(_worker, t)
            in_flight_tasks[fut] = t
        print(f"    [SUBMIT] Initial window ready. Entering drain loop...")

        while in_flight_tasks:
            # Block until at least one future completes
            done, _ = wait(in_flight_tasks.keys(), return_when=FIRST_COMPLETED)

            for fut in done:
                res  = fut.result()
                task = in_flight_tasks.pop(fut)
                completed += 1

                # Submit next task immediately to keep workers busy
                t = next(task_iter, None)
                if t is not None:
                    new_fut = executor.submit(_worker, t)
                    in_flight_tasks[new_fut] = t

                if completed % 500 == 0 or completed == total:
                    print(f"    Progress: {completed}/{total} grid points done "
                        f"({skipped + completed}/{N1*N2*N3*N4} total).")

                if res is None:
                    continue

                i, j, l, m = res['i'], res['j'], res['l'], res['m']

                if 'Logan_Vt' in res['true_values']:
                    mm_true[i, j, l, m, 0] = res['true_values']['Logan_Vt']
                if 'Patlak_Ki' in res['true_values']:
                    mm_true[i, j, l, m, 1] = res['true_values']['Patlak_Ki']

                mm_params[i, j, l, m, 0] = task['sim_params']['K1']
                mm_params[i, j, l, m, 1] = task['sim_params']['k2']
                mm_params[i, j, l, m, 2] = task['sim_params']['k3']
                mm_params[i, j, l, m, 3] = task['sim_params']['k4']

                for nl, ts_dict in res['noise_tstart_results'].items():
                    for ts, mt_dict in ts_dict.items():
                        mm_tac[(nl, ts)][i, j, l, m, 0, :] = res['noise_tstart_tac'][nl][ts]['time']
                        mm_tac[(nl, ts)][i, j, l, m, 1, :] = res['noise_tstart_tac'][nl][ts]['tac']

                        for model_type, metrics in mt_dict.items():
                            mm[(nl, ts, model_type)][i, j, l, m, 0]  = metrics['slope']
                            mm[(nl, ts, model_type)][i, j, l, m, 1]  = metrics['ic']
                            mm[(nl, ts, model_type)][i, j, l, m, 2]  = metrics['error']
                            mm[(nl, ts, model_type)][i, j, l, m, 3:] = metrics['model']

                mm_sen[i, j, l, m] = True

    # Flush all memmaps explicitly before returning
    mm_true.flush()
    mm_sen.flush()

    for arr in mm_tac.values():
        arr.flush()

    for arr in mm.values():
        arr.flush()

# ─────────────────────────────────────────────────────────────────────────────
# ENTRY POINTS
# ─────────────────────────────────────────────────────────────────────────────

def run(K1, k2, k3, k4, vb, methods, datapoints, organ,
        noise_values, tstart_values, framebased_bif, interpolated_bif,
        base):

    npys = Path(base) / "npys"
    temp = Path(base) / "temp"
    orgDir = Path(npys) / organ
    
    npys.mkdir(parents=True, exist_ok=True)
    temp.mkdir(parents=True, exist_ok=True)
    orgDir.mkdir(parents=True, exist_ok=True)

    framebased_bif_min = framebased_bif.replace('.txt','_minutes.txt')
    _convert_sec2min(framebased_bif,framebased_bif_min)
    simulate_and_fit(
        K1_list          = K1,
        k2_list          = k2,
        k3_list          = k3,
        k4_list          = k4,
        vb               = vb,
        methods          = methods,
        framebased_bif   = framebased_bif_min,
        interpolated_bif = interpolated_bif,
        noise_list       = list(noise_values),
        datapoints       = datapoints,
        tstart_seconds   = [int(ts) for ts in tstart_values],
        n_workers        = 9,
        orgDir           = str(orgDir),
        temp_dir         = str(temp),
    )

# ─────────────────────────────────────────────────────────────────────────────
# DEFINE TISSUE VALUES
# ─────────────────────────────────────────────────────────────────────────────

def define_tissue_values(tissue):
    if tissue   == 'tumour':       tissue_values = [[11.9, 3],    [0.112, 0.05],  [0.304, 0.1],   [0.140, 0.07],  [0.005, 0.005]]
    elif tissue == 'liver':        tissue_values = [[0.1,0.2],   [0.667, 0.223],   [0.695, 0.175],   [0.04, 0.005], [0.02, 0.002]]
    elif tissue == 'lungs':        tissue_values = [[14.5, 8.7],  [0.233, 0.447], [1.536, 2.305], [0.012, 0.005], [0.005, 0.004]]
    elif tissue == 'kidneys':      tissue_values = [[14.6, 8],    [0.793, 0.277], [0.672, 0.324], [0.008, 0.006], [0.004, 0.002]]
    elif tissue == 'spleen':       tissue_values = [[4.8, 4.2],   [1.535, 0.546], [2.697, 0.496], [0.007, 0.003], [0.004, 0.002]]
    elif tissue == 'bones':        tissue_values = [[2, 1],       [0.108, 0.056], [0.642, 0.195], [0.025, 0.004], [0.01, 0.01]]
    elif tissue == 'grey_matter':  tissue_values = [[4.2, 2],     [0.130, 0.026], [0.159, 0.066], [0.064, 0.031], [0.02, 0.02]]
    elif tissue == 'white_matter': tissue_values = [[2.1, 1],     [0.064, 0.015], [0.090, 0.031], [0.026, 0.013], [0.01, 0.01]]
    elif tissue == 'Irreversible': tissue_values = [[2.1, 1],     [0.064, 0.015], [0.090, 0.031], [0.026, 0.013], [0.00, 0.00]]
    
    try:
        return tissue_values
    except NameError:
        print('Tissue not recognised. Please choose from: tumour, liver, lungs, '
              'kidneys, spleen, bones, grey_matter or white_matter.')

