"""
-----------------------------------------------------------------------
Script for Unsteady Aerodynamics using the Doublet-Lattice Method (DLM)

Author: Higor Luis Silva
-----------------------------------------------------------------------
"""
import sys
sys.path.append('..')
import numpy as np
from pyaero.wing import WingParameters
from pyaero.discretize import discretize
from pyaero.calcD import DLMParameters, calcD

def test_DLM():
    paramDLM = DLMParameters()
    wing = WingParameters()

    totalboxes, box = discretize(wing, paramDLM)
    k = 1
    D = calcD(wing, box, totalboxes, k, paramDLM=paramDLM)
    AIC = np.linalg.inv(D)
    return 0
