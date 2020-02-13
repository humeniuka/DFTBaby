# -*- coding: utf-8 -*-
"""

References
----------
[1]  Jorge Nocedal and Stephen Wright, 
     "Numerical Optimization", 
     Springer

"""
import numpy as np

def log_barrier(C, A):
    """
    The log-barrier prevents the parameter vector from leaving the feasible area, where the constraints
    are fulfilled by adding a repulsive barrier B(x) to the objective function f(x):

      B(x) = - sum_i log(C_i(x))

    The gradient of the barrier becomes

      dB/dxj = - sum_i 1/C_i(x) dC_i/dxj = - sum_i 1/C_i A_ij

    The auxiliary function, that is minimized instead of f(x), is now

      f(x) + nu*B(x)

    where ``nu`` is a small adjustable number.

    Parameters
    ----------
    C: vector with values of constraints C_i(x), i=1,...,m
    A: (m x n) matrix, where A[i,j] contains the derivative of the i-th constraint w/r/t the j-th variable

    Returns
    -------
    B: value of barrier function
    dB: gradient of barrier function
    """
    m = len(C)   # number of constraints
    B = 0.0
    for i in range(0, m):
        B -= np.log(C[i])
    # gradient
    dB = - np.dot(1.0/C, A)

    return B, dB

