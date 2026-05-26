######################################################################################
##################################### Packages  ######################################
######################################################################################
import os
import shutil
from datetime import datetime
from subprocess import Popen
import numpy as np

import petkinetics.bin.py.Patlak.Patlak_fits as patlak
import petkinetics.bin.py.dLogan.dLogan_fits as logan
import petkinetics.bin.py.getPath as get_path
from petkinetics.prototype_bridge import bridge

# Global constants:

kiStart_glob = 0.005
vdStart_glob = 0.1
kiBounds_glob = [0.0001, 0.5]
vdBounds_glob = [0.0001, 2.0]

K1start_glob = 0.1
K2start_glob = 0.1
K1bounds_glob = [0.01,2]
K2bounds_glob = [0.01,20]

alphastart_glob = 1.0
betarstart_glob = -10.0
alphabounds_glob = [0.0, 100.0]
betabounds_glob = [-1000, -0.01]

######################################################################################
############################ Region Based Logan Fitting  #############################
######################################################################################
'''
tacFile should be formatted as: [frame time (integer seconds)] [tissue value (Bq/mL)] 
and should include only the steady-state frames.

bifFile should be formatted as: [time (integer seconds)] [blood value (Bq/mL)]
and should include the entire input function from 0s until the end of the scan,
in 1s resolution.
'''

def Patlak_MultiStar(
        tacFile : str, 
        bifFile : str, 
        numStars : int, 
        kiBounds : int = None, 
        vdBounds : np.ndarray = None, 
        figFlag : bool =True):
    """
    Fits Ki and Vd using Patlak multistar prototype

    Parameters
    ----------
    tacFile : str
        Path to textfile with time-activity curve.
    bifFile : str
        Path to textfile with blood input-function.
    numStars : int
        Number of star-values
    kiBounds : int  
        Array of ki bounds, e.g. [0.0001, 2]
    vdBounds : np.ndarray 
        Array of vd bounds, e.g. [0.0001, 4]
    figFlag : bool
        Show figure of the fit. 
    
    Returns
    -------
    ki : float
        Fitted Ki 
    vd : float
        fitted vd
    error : float
        Sum of squared errors 
    model : np.ndarray 
        The fitted model
    frameTime : np.ndarray 
        The fitted frame times
    tac : np.ndarray
        The fitted activity values
    """
    
    if not kiBounds and not vdBounds:
        kiBounds, vdBounds = kiBounds_glob, vdBounds_glob

    return patlak.Patlak_MultiStar(tacFile, bifFile, numStars, kiBounds, vdBounds, figFlag=figFlag)

def Patlak_tStar(tacFile : str, 
        bifFile : str, 
        kiStart : float = None, 
        vdStart : float = None,
        kiBounds : int = None, 
        vdBounds : np.ndarray = None, 
        figFlag : bool =True,
        linear : bool = True):
    """
    Fits Ki and Vd using Patlak tStar prototype

    Parameters
    ----------
    tacFile : str
        Path to textfile with time-activity curve.
    bifFile : str
        Path to textfile with blood input-function.
    kiStart : float
        Initial ki guess
    vdStart : float
        Initial vd guess
    kiBounds : np.ndarray   
        Array of ki bounds, e.g. [0.0001, 2]
    vdBounds : np.ndarray 
        Array of vd bounds, e.g. [0.0001, 4]
    figFlag : bool
        Show figure of the fit. 
    
    Returns
    -------
    ki : float
        Fitted Ki 
    vd : float
        fitted vd
    error : float
        Sum of squared errors 
    model : np.ndarray 
        The fitted model
    frameTime : np.ndarray 
        The fitted frame times
    tac : np.ndarray
        The fitted activity values
    """
    if not kiStart and not vdStart and not kiBounds and not vdBounds:
        kiStart, vdStart, kiBounds, vdBounds = kiStart_glob, vdStart_glob, kiBounds_glob, vdBounds_glob

    return patlak.Patlak_tStar(tacFile, bifFile, kiStart, vdStart, kiBounds, vdBounds, figFlag=figFlag, linear=linear)

def Patlak_tZero(tacFile : str, 
        bifFile : str, 
        kiStart : float = None, 
        vdStart : float = None,
        kiBounds : int = None, 
        vdBounds : np.ndarray = None, 
        figFlag : bool =True):
    """
    Fits Ki and Vd using Patlak tZero prototype

    Parameters
    ----------
    tacFile : str
        Path to textfile with time-activity curve.
    bifFile : str
        Path to textfile with blood input-function.
    kiStart : float
        Initial ki guess
    vdStart : float
        Initial vd guess
    kiBounds : np.ndarray   
        Array of ki bounds, e.g. [0.0001, 2]
    vdBounds : np.ndarray 
        Array of vd bounds, e.g. [0.0001, 4]
    figFlag : bool
        Show figure of the fit. 
    
    Returns
    -------
    ki : float
        Fitted Ki 
    vd : float
        fitted vd
    error : float
        Sum of squared errors 
    model : np.ndarray 
        The fitted model
    frameTime : np.ndarray 
        The fitted frame times
    tac : np.ndarray
        The fitted activity values
    """
    if not kiStart and not vdStart and not kiBounds and not vdBounds:
        kiStart, vdStart, kiBounds, vdBounds = kiStart_glob, vdStart_glob, kiBounds_glob, vdBounds_glob

    return patlak.Patlak_tZero(tacFile, bifFile, kiStart, vdStart, kiBounds, vdBounds, figFlag=figFlag)

def dLoganRB_tStar(tacFile : str, 
                   bifFile : str, 
                   slp_start : float = None, 
                   ic_start : float = None, 
                   slp_blounds : np.ndarray = None, 
                   ic_bounds : np.ndarray = None, 
                   k1k2Flag : bool = True, 
                   figFlag : bool =False):
    """
    Fits Vt and Beta using Logan tStar prototype

    Parameters
    ----------
    tacFile : str
        Path to textfile with time-activity curve.
    bifFile : str
        Path to textfile with blood input-function.
    slp_start: float
        Initial guess of slope (Vt)
    ic_start : float
    slp_blounds : np.ndarray 
        Array of vt bounds, e.g. [0.0, 20]
    ic_bounds : np.ndarray 
        Array of beta bounds, e.g. [-1000, -0]
    k1k2Flag: bool
        Defines if model is based on 1-tissue compartment model (True) or 2-tissue compartment model (False)
    figFlag : bool
        Show figure of the fit. 
    
    Returns
    -------
    vt : float
        Fitted Vt 
    beta : float
        fitted beta
    error : float
        Sum of squared errors 
    model : np.ndarray 
        The fitted model
    frameTime : np.ndarray 
        The fitted frame times
    tac : np.ndarray
        The fitted activity values
    """
    if not slp_start and not ic_start and not slp_blounds and not ic_bounds:
        if k1k2Flag:
            slp_start, ic_start, slp_blounds, ic_bounds = K1start_glob, K2start_glob, K1bounds_glob, K2bounds_glob
        if not k1k2Flag:
            slp_start, ic_start, slp_blounds, ic_bounds = alphastart_glob, betarstart_glob, alphabounds_glob, betabounds_glob

    return logan.dLoganRB_tStar(tacFile, bifFile, slp_start, ic_start, slp_blounds, ic_bounds, k1k2Flag = k1k2Flag, figFlag=figFlag)

def dLoganRB_tZero(tacFile : str, 
                   bifFile : str, 
                   slp_start : float = None, 
                   ic_start : float = None, 
                   slp_blounds : np.ndarray = None, 
                   ic_bounds : np.ndarray = None, 
                   k1k2Flag : bool = True, 
                   figFlag : bool =False):
    """
    Fits Vt and Beta using Logan tStar prototype

    Parameters
    ----------
    tacFile : str
        Path to textfile with time-activity curve.
    bifFile : str
        Path to textfile with blood input-function.
    slp_start: float
        Initial guess of slope (Vt)
    ic_start : float
    slp_blounds : np.ndarray 
        Array of vt bounds, e.g. [0.0, 20]
    ic_bounds : np.ndarray 
        Array of beta bounds, e.g. [-1000, -0]
    k1k2Flag: bool
        Defines if model is based on 1-tissue compartment model (True) or 2-tissue compartment model (False)
    figFlag : bool
        Show figure of the fit. 
    
    Returns
    -------
    vt : float
        Fitted Vt 
    beta : float
        fitted beta
    error : float
        Sum of squared errors 
    model : np.ndarray 
        The fitted model
    frameTime : np.ndarray 
        The fitted frame times
    tac : np.ndarray
        The fitted activity values
    """
    if not slp_start and not ic_start and not slp_blounds and not ic_bounds:
        if k1k2Flag:
            slp_start, ic_start, slp_blounds, ic_bounds = K1start_glob, K2start_glob, K1bounds_glob, K2bounds_glob
        if not k1k2Flag:
            slp_start, ic_start, slp_blounds, ic_bounds = alphastart_glob, betarstart_glob, alphabounds_glob, betabounds_glob

    return logan.dLoganRB_tZero(tacFile, bifFile, slp_start, ic_start, slp_blounds, ic_bounds, k1k2Flag = k1k2Flag, figFlag=figFlag)


######################################################################################
################################ Voxel Based Fitting  ################################
######################################################################################

## Internal helpers:

def _create_dir(src):
    """
    Creates directory if it does not already exists.
    """
    if not os.path.isdir(src):
        os.makedirs(src)

def _prep_parametric_folder(src):
    
    """
    Organises parametric directory if necessary.
    """

    approved_dirs = ["bif.txt", "dicoms", "sliceTimes", "volumes", "Logan", "Patlak","prev", "SUV", "frameTimes.txt","bifSS.txt"]
    lst = os.listdir(src)
    fold = []
    fold = [obj for obj in lst if not [app_dir for app_dir in approved_dirs if app_dir in obj]]

    if fold:
        dst = os.path.join(src,'prev')
        _create_dir(dst)
        x = (datetime.now()).strftime("%d-%b-%Y-%H.%M.%S")

        for dir in fold: 
            src_dir = os.path.join(src,dir)

            if os.path.isfile(src_dir):

                suffix = (dir.split('.'))[-1]
                base_name = ''

                for part in (dir.split('.')[0:-1]):
                    base_name = base_name+part+'_'

                dst_dir = os.path.join(dst,base_name+x+'.'+suffix)
            if os.path.isdir(src_dir):
                dst_dir = rf'{os.path.join(dst,dir)}-{x}'

            shutil.move(src_dir,dst_dir)

def _run_exe(modelExe, src, initVals, slpLim, intcptLim, itr, flag, opt):
    """
    Runs executable for parametric image and generate log.
    """
    log = os.path.join(src, "log.txt")       # Path for patient file
    ITF = os.path.join(src, "ITF.txt")       # Path for patient file

    with open(ITF, 'w') as f:
        f.write(src+'\n'+initVals+'\n'+ slpLim +'\n'+ intcptLim +'\n'+ itr +'\n'+ flag)  
    command = modelExe + " " + ITF + " " + opt

    with open(log, 'w') as f:
        f.write(f"Timestamp: {datetime.now()}\n\nExecutable: {modelExe}\nInitial values: {initVals}\nSlope Limits: {slpLim}\nIntercept Limits: {intcptLim}\nIterations: {itr}\n\nCommand: {command}")  

    process = Popen(command)
    process.communicate()   

def _export_err2dcm(src):
    """
    Exports error from .raw image to dicom
    """
    for folders in os.listdir(src): 
        current = os.path.join(src, folders)
        if 'Error' in current and '.raw' in current and os.path.isfile(current):
            raw = current
        if 'Beta_' in current and os.path.isdir(current):
            dicom = current
        if 'Vd_' in current and os.path.isdir(current):
            dicom = current

        if raw and dicom: continue
    
    bridge.exportErr2dcm(raw,dicom)

def _organise_parametric_folder(src, model, type, cm,  alphaStart, betaStart):

    """
    Organises new parametric folders and files.
    """

    approved_dirs = ["bif.txt", "dicoms", "sliceTimes", "volumes", "Logan", "Patlak","prev", "SUV", "frameTimes.txt","bifSS.txt"]
    lst = os.listdir(src)
    fold = []
    fold = [obj for obj in lst if not [app_dir for app_dir in approved_dirs if app_dir == obj]]

    if fold:
        dst = os.path.join(src, model); _create_dir(dst)
        dst_t = os.path.join(dst, type); _create_dir(dst_t)

        x = (datetime.now()).strftime("%d-%b-%Y")

        dst_folder = os.path.join(dst_t, f'{cm}_{alphaStart}_{betaStart}_{x}'); _create_dir(dst_folder)

        for dir in fold: 
            src_dir = os.path.join(src,dir)

            if os.path.isfile(src_dir):

                suffix = (dir.split('.'))[-1]
                base_name = ''

                for part in (dir.split('.')[0:-1]):
                    base_name = base_name+part+'_'

                dst_dir = os.path.join(dst_folder,base_name+'.'+suffix)

            if os.path.isdir(src_dir):
                dst_dir = rf'{os.path.join(dst_folder,dir)}'

            shutil.move(src_dir,dst_dir)

def _run_parametric_model(modelExe, src, alphaStart, betaStart, slpLim, intcptLim, itr, model, type, cm, errFlg=True, lm = False):

    """ 
    Runs parametric model 
    """
    _prep_parametric_folder(src)

    initVals = f'{str(alphaStart)} {str(betaStart)}'
    bif = os.path.join(src, 'bif.txt') 
    dicoms = os.path.join(src, 'dicoms') 
    
    if errFlg == True: flag = '1'
    elif errFlg == False: flag = '0'

    if lm: opt = "-d -lm"
    else:  opt = "-d"
    
    if os.path.isfile(bif) and os.path.isdir(dicoms):
        
        _run_exe(modelExe, src, initVals, slpLim, intcptLim, itr, flag, opt)

        if errFlg: 
            _export_err2dcm(src)

        _organise_parametric_folder(src, model, type, cm,  alphaStart, betaStart)
    else: 
        print("Blood input function or dicom folder doesn't exist.")
    

## Run parametric model

def PatlakVB_Tzero_CBM(inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr = 30, errFlg=True, lm = False):
    modelExe = get_path.getVB_Patlak_Tzero_CBM()
    model       = 'Patlak'
    type        = 'tZero'
    cm          = 'T2_CBM'

    _run_parametric_model(modelExe, inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr, model, type, cm, errFlg)

def PatlakVB_Tstar_SBP(inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr = 30, errFlg=True, lm = False):

    if lm: modelExe = get_path.getVB_Patlak_Tstar_SBP_LM()
    else: modelExe = get_path.getVB_Patlak_Tstar_SBP_GD()

    model       = 'Patlak'
    type        = 'tStar'
    cm          = 'T2_SBP'

    _run_parametric_model(modelExe, inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr, model, type, cm, errFlg)

def dLoganVB_Tzero_CBM(inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr = '30', errFlg=True, lm=False):
    modelExe    = get_path.getVB_Logan_Tzero_CBM()
    model       = 'Logan'
    type        = 'tZero'
    cm          = 'T2_CBM'

    _run_parametric_model(modelExe, inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr, model, type, cm, errFlg, lm=lm)

def dLoganVB_Tzero_1T_CBM(inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr='30', errFlg=True, lm=False):
    modelExe    = get_path.getVB_Logan_Tzero_T1_CBM()
    model       = 'Logan'
    type        = 'tZero'
    cm          = 'T1_CBM'

    _run_parametric_model(modelExe, inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr, model, type, cm, errFlg, lm=lm)

def dLoganVB_Tzero_1T_SBP(inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr='30', errFlg=True):
    modelExe = get_path.getVB_Logan_Tzero_T1_SBP()
    model       = 'Logan'
    type        = 'tZero'
    cm          = 'T1_SBP'

    _run_parametric_model(modelExe, inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr, model, type, cm, errFlg)

def dLoganVB_Tstar_1T_SBP(inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr='30', errFlg=True):
    modelExe = get_path.getVB_Logan_Tstar_T1_SBP()
    model       = 'Logan'
    type        = 'tStar'
    cm          = 'T1_SBP'

    _run_parametric_model(modelExe, inputDir, alphaStart, betaStart, alphaBounds,betaBounds, itr, model, type, cm, errFlg)
