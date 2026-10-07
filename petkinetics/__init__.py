"""
petkinetics — PET Kinetic Modelling Toolkit
============================================
A companion package to ``medimkit`` for kinetic modelling of dynamic PET data.
Covers linearised graphical analysis, model comparison, TAC simulation, and
bridging to proprietary prototype models.

Submodules
----------
model_fitting       : Fit Logan and Patlak graphical models to TAC data
model_selection     : Compare models via goodness-of-fit and statistical criteria
tac_simulation      : Simulate time-activity curves from known model parameters
prototype_bridge    : Call proprietary executables and DLL-based models
"""

__version__ = "0.1.0"
__author__ = "Josefine Rosenskjold Madsen"

from petkinetics import (
    model_fitting,
    model_selection,
    simulation,
    prototype_bridge,
    input_functions,

)

__all__ = [
    "model_fitting",
    "model_selection",
    "tac_simulation",
    "prototype_bridge",
    "input_functions",
]