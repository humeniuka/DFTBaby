#!/usr/bin/env python

"""
In order to use user-defined CVs, modify the following functions and include this file in your \
working directory.
"""

import numpy as np

def customcv(symb, coords, prm):
    """
    Return value of user-defined CV for given symbols, coordinates and prm
    
    Parameters:
    -----------
    symb: list of strings
        atomic symbols of the molecule
    coords: numpy array
        coordinates of the molecule (shape: (Nat, 3))
    prm: dict
        parameters given in the CV section of meta-config.py
    
    Returns:
    --------
    custom: number
        CV value
    
    Note:
    -----
    Do not modify function name or arguments.
    
    """
    ### modify the following code
    custom = 0.0
    ### end modify
    return custom

def dcustomcv(symb, coords, prm):
    """
    Return derivative of user_defined CV wrt coordinates for given symbols, coordinates and prm
    
    Parameters:
    -----------
    symb: list of strings
        atomic symbols of the molecule
    coords: numpy array
        coordinates of the molecule (shape: (Nat, 3))
    prm: dict
        parameters given in the CV section of meta-config.py
    
    Returns:
    --------
    dcustom: numpy array
        derivative of CV wrt coordinates (shape: (Nat, 3))
    
    Note:
    -----
    Do not modify function name or arguments.
    
    """
    ### modify the following code
    dcustom = np.zeros((len(symb), 3))
    ### end modify
    return dcustom
