import numpy as np
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap


def calculate_AIC(model1_sse: np.ndarray, model2_sse: np.ndarray, n: int, k1: int, k2: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute voxelwise Akaike Information Criterion (AIC) for two competing models.

    AIC balances goodness-of-fit against model complexity, penalising models
    with more free parameters. It is calculated as:

        AIC = n * log(SSE / n) + 2k

    where SSE is the sum of squared errors, n is the number of data points
    (time frames), and k is the number of free parameters in the model.

    Parameters
    ----------
    model1_sse : np.ndarray
        Voxelwise sum of squared errors for model 1. Shape should match
        model2_sse.
    model2_sse : np.ndarray
        Voxelwise sum of squared errors for model 2.
    n : int
        Number of time frames used in model fitting.
    k1 : int
        Number of free parameters in model 1. 
    k2 : int
        Number of free parameters in model 2. 
    
    Returns
    -------
    model1_aic : np.ndarray
        Voxelwise AIC values for model 1.
    model2_aic : np.ndarray
        Voxelwise AIC values for model 2.

    Notes
    -----
    This implementation uses the standard AIC formulation and assumes
    normally distributed residuals. For small sample sizes (n/k < 40),
    consider using the corrected form AICc = AIC + 2k(k+1)/(n-k-1).
    """
    model1_aic = n * np.log(model1_sse / n) + 2 * k1
    model2_aic = n * np.log(model2_sse / n) + 2 * k2
    return model1_aic, model2_aic

def calculate_AICc(model1_sse: np.ndarray, model2_sse: np.ndarray, n: int, k1: int, k2: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute voxelwise corrected Akaike Information Criterion (AICc) for two
    competing models.

    AICc extends AIC with a correction term for small sample sizes:

        AICc = n * log(SSE / n) + 2k + 2k(k+1) / (n - k - 1)

    The correction term becomes negligible as n grows large relative to k,
    and AICc converges to AIC. As a rule of thumb, AICc is preferred when
    n/k < 40; for larger ratios the two criteria give effectively identical
    results.

    Parameters
    ----------
    model1_sse : np.ndarray
        Voxelwise sum of squared errors for model 1. Shape must match
        model2_sse.
    model2_sse : np.ndarray
        Voxelwise sum of squared errors for model 2.
    n : int
        Number of time frames used in model fitting.
    k1 : int
        Number of free parameters in model 1. 
    k2 : int
        Number of free parameters in model 2. 

    Returns
    -------
    model1_aicc : np.ndarray
        Voxelwise AICc values for model 1.
    model2_aicc : np.ndarray
        Voxelwise AICc values for model 2.

    Raises
    ------
    ValueError
        If n <= k + 1, making the correction term undefined (division by zero).

    """
    if n <= k1 + 1 or n <= k2 + 1:
        raise ValueError(
            f"n must be greater than k + 1 for AICc to be defined, "
            f"got n={n}, k1={k1} and k2={k2}. Use calculate_AIC() instead, or "
            f"increase the number of time frames."
        )

    c1 = (2 * k1 * (k1 + 1)) / (n - k1 - 1)
    c2 = (2 * k2 * (k2 + 1)) / (n - k2 - 1)
    model1_aicc = n * np.log(model1_sse / n) + 2 * k1 + c1
    model2_aicc = n * np.log(model2_sse / n) + 2 * k2 + c2
    return model1_aicc, model2_aicc

def calculate_delta_map(model1_sse: np.ndarray, model2_sse: np.ndarray, n: int, k: int, corrected_form = False) -> np.ndarray:
    """
    Compute a voxelwise ΔAIC map between two competing kinetic models.

    ΔAIC = AIC(model1) − AIC(model2)

    Negative values indicate model 1 is preferred; positive values indicate
    model 2 is preferred. The magnitude reflects the strength of evidence:
    |ΔAIC| > 2 is generally considered meaningful, and |ΔAIC| > 10 is strong.

    Parameters
    ----------
    model1_sse : np.ndarray
        Voxelwise sum of squared errors for model 1.
    model2_sse : np.ndarray
        Voxelwise sum of squared errors for model 2. Shape must match
        model1_sse.
    n : int
        Number of time frames used in model fitting.
    k : int
        Number of free parameters. Must be the same for both models.
    corrected_form : bool
        Defining if AIC should be in corrected form. Default is False.


    Returns
    -------
    np.ndarray
        Voxelwise ΔAIC map (model1 − model2). Negative values favour model 1;
        positive values favour model 2.

        
    """
    model1_aic, model2_aic = calculate_AIC(model1_sse, model2_sse, n, k)
    return np.subtract(model1_aic, model2_aic)

def createErrorMap(model1_sse: np.ndarray, model2_sse: np.ndarray) -> np.ndarray:
    """
    Compute a voxelwise SSE difference map between two models.

    Returns the elementwise difference in sum of squared errors:

        error_map = SSE(model1) − SSE(model2)

    Negative values indicate model 1 has lower residual error at that voxel;
    positive values indicate model 2 fits better. Unlike ΔAIC, this map does
    not penalise for model complexity and should be interpreted alongside a
    model selection criterion such as ΔAIC.

    Parameters
    ----------
    model1_sse : np.ndarray
        Voxelwise sum of squared errors for model 1.
    model2_sse : np.ndarray
        Voxelwise sum of squared errors for model 2. Shape must match
        model1_sse.

    Returns
    -------
    np.ndarray
        Voxelwise SSE difference map (model1 − model2).

    Notes
    -----
    This map is sensitive to absolute fit quality but does not account for
    model complexity. It is most useful as a diagnostic alongside ΔAIC —
    for example, to identify regions where one model fits substantially
    better regardless of complexity penalty.
    """
    return np.subtract(model1_sse, model2_sse)

def delta_cmap(scale = 7, threshold = 2, bad=False):
    """
    Create a custom colormap for delta maps with a white center.
    
    Parameters:
    vmin (float): Minimum value for the colormap.
    vmax (float): Maximum value for the colormap.
    mid (float): Center value for the colormap.
    threshold (float): Threshold around zero to be white.
    
    Returns:
    custom_cmap: A LinearSegmentedColormap object.
    norm: A TwoSlopeNorm object for normalization.
    """
    part = threshold/(scale*2)

    colormap = LinearSegmentedColormap.from_list('custom',
                                                [(0,'#0000FF'),
                                                (0.5-part,'#ffffff'),
                                                (0.5+part,'#ffffff'),
                                                (1,'#ff0000')],256)

    norm = colors.TwoSlopeNorm(vmin=-scale, vcenter=0, vmax=scale)

    if bad:
        colormap.set_bad(color='g')  # Set color for NaN values
        

    return colormap, norm