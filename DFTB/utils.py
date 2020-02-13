"""
odd bits and pieces of code
"""
from numpy import zeros, argsort
import numpy as np
import inspect

# for python >= 2.7 this function is added in itertools
def combinations_with_replacement(iterable, r):
    # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
    pool = tuple(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        yield tuple(pool[i] for i in indices)

# for scipy >= v0.11 this function is added to scipy.linalg
# This implementation is stolen from scipy v0.11 linalg/special_matrices.py
def block_diag(*arrs):
    """
    Create a block diagonal matrix from provided arrays.

    Given the inputs `A`, `B` and `C`, the output will have these
    arrays arranged on the diagonal::

        [[A, 0, 0],
         [0, B, 0],
         [0, 0, C]]

    Parameters
    ----------
    A, B, C, ... : array_like, up to 2-D
        Input arrays.  A 1-D array or array_like sequence of length `n`is
        treated as a 2-D array with shape ``(1,n)``.

    Returns
    -------
    D : ndarray
        Array with `A`, `B`, `C`, ... on the diagonal.  `D` has the
        same dtype as `A`.

    Notes
    -----
    If all the input arrays are square, the output is known as a
    block diagonal matrix.

    Examples
    --------
    >>> from scipy.linalg import block_diag
    >>> A = [[1, 0],
    ...      [0, 1]]
    >>> B = [[3, 4, 5],
    ...      [6, 7, 8]]
    >>> C = [[7]]
    >>> block_diag(A, B, C)
    [[1 0 0 0 0 0]
     [0 1 0 0 0 0]
     [0 0 3 4 5 0]
     [0 0 6 7 8 0]
     [0 0 0 0 0 7]]
    >>> block_diag(1.0, [2, 3], [[4, 5], [6, 7]])
    array([[ 1.,  0.,  0.,  0.,  0.],
           [ 0.,  2.,  3.,  0.,  0.],
           [ 0.,  0.,  0.,  4.,  5.],
           [ 0.,  0.,  0.,  6.,  7.]])

    """
    if arrs == ():
        arrs = ([],)
    arrs = [np.atleast_2d(a) for a in arrs]

    bad_args = [k for k in range(len(arrs)) if arrs[k].ndim > 2]
    if bad_args:
        raise ValueError("arguments in the following positions have dimension "
                            "greater than 2: %s" % bad_args)

    shapes = np.array([a.shape for a in arrs])
    out = np.zeros(np.sum(shapes, axis=0), dtype=arrs[0].dtype)

    r, c = 0, 0
    for i, (rr, cc) in enumerate(shapes):
        out[r:r + rr, c:c + cc] = arrs[i]
        r += rr
        c += cc
    return out

def argsort_2d(arr):
    """
    for a 2-dimensional array return two lists of indeces
    row_sort and col_sort such that arr[row_sort[i],col_sort[i]] is the
    i-th lowest element in arr.

    Parameters:
    ===========
    arr: 2D numpy array
    
    Returns:
    ========
    row_sort: 1D numpy array with row indeces
    col_sort= 1D numpy array with column indeces
    """
    assert len(arr.shape) == 2
    nrows,ncols = arr.shape
    sort_indx = argsort(arr, axis=None)
    row_sort = sort_indx / ncols
    col_sort = sort_indx % ncols
    
    return row_sort, col_sort

def annotated_matrix(M, row_labels, col_labels, format="%.4f", colwidth=10, block_len=10):
    """
    format a matrix such that columns and rows are labeled. If rows are too
    long they are split into blocks.

    Parameters:
    ===========
    M: 2d numpy array with shape mxn
    row_labels: list of m labels for rows
        a label "-" creates a horizontal separation line
    col_labels: list of n labels for columns

    Returns:
    ========
    string with formatted matrix
    """
    import string
    from math import ceil
    m,n = M.shape

    nblocks = int(ceil(n/float(block_len)))
    txt = ""
    for b in range(0, nblocks):
        txt += " "*(colwidth+1) + "|"
        for col in range(b*block_len, min(n,(b+1)*block_len)):
            txt += string.center(col_labels[col], colwidth+1) + "|"
        txt += "\n" + "-"*(colwidth+2)*(min(block_len,n)+1) + "\n"

        nr_sep = 0 # count the separation lines to keep track
                   # of how many data lines were printed
        for row in range(0,len(row_labels)):
            if row_labels[row] == "-":
                txt += "-"*(min(block_len,n)+1)*(colwidth+2)
                nr_sep += 1
            else:
                txt += string.center(row_labels[row], colwidth+1) + "|"
                for col in range(b*block_len, min(n,(b+1)*block_len)):
                    txt += string.center(format % M[row-nr_sep,col], colwidth+1) + "|"
            txt += "\n"
        txt += "\n"
    return txt


class dotdic(dict):
    """
    overload dictionary to allow accessing data by .-notation
    e.g.
    >> d = dotdic()
    >> d["bla"] = 1
    >> print d.bla
    """
    def __getattr__(self, name):
        return self[name]
    def __setattr__(self, name, value):
        self[name] = value


def numerical_gradient(f,x0,h=1.0e-5):
    """
    compute gradient of f at x0 by numerical differentiation

    Parameters:
    ===========
    f: scalar function of a 1D vector x
    x0: 1D numpy array

    Returns:
    ========
    1D numpy array with difference quotient 
       df/dx_i = (f(x0 + h*ei) - f(x0))/h 
    """
    n = len(x0)
    f0 = f(x0)
    dfdx = zeros(n)
    for i in range(0, n):
        print "numerical gradient: %d of %d" % (i,n)
        ei = zeros(n)
        ei[i] = 1.0 # unit vector
#        # forward gradient
#        x_hi = x0 + h*ei
#        dfdx[i] = (f(x_hi) - f0)/h
        # symmetric gradient
        x_mhi = x0 - h*ei
        x_phi = x0 + h*ei
        dfdx[i] = (f(x_phi) - f(x_mhi))/(2.0*h)
    return dfdx

def numerical_hessian(f,x0,h=1.0e-8):
    """
    compute matrix of second derivatives of f at f0 by finite differences

    Parameters:
    ===========
    f: scalar function of a 1D vector x
    x0: 1D numpy array

    Returns:
    ========
    2D numpy array with difference quotient 
       H_ij = d^2f/(dx_i dx_j)(x0)

    Warning: This is probably not the most efficient way to compute second derivatives,
    since f(x) is evaluated several times at the same x.
    """
    n = len(x0)
    hessian = zeros((n,n))
    f0 = f(x0)
    for i in range(0, n):
        ei = zeros(n)
        ei[i] = 1.0
        for j in range(i, n):
            if i != j:
                ej = zeros(n)
                ej[j] = 1.0
                x_pipj = x0 + h*(ei + ej)
                x_pimj = x0 + h*(ei - ej)
                x_mipj = x0 + h*(-ei + ej)
                x_mimj = x0 - h*(ei + ej)

                hessian[i,j] = ( f(x_pipj) - f(x_pimj) - f(x_mipj) + f(x_mimj) )/(4*h*h)
                hessian[j,i] = hessian[i,j]
        # i == j
        x_pi = x0 + h*ei
        x_mi = x0 - h*ei
        hessian[i,i] = ( f(x_pi) - 2*f0 + f(x_mi) )/(h*h)
    return hessian

def numerical_hessian_G(grad,x0,h=1.0e-8):
    """
    compute hessian by numerical differentiation of gradients

    Parameters:
    ===========
    grad: grad(x) computes the gradient at position x
    """
    n = len(x0)
    hessian = zeros((n,n))
    g0 = grad(x0)
    for i in range(0, n):
        ei = zeros(n)
        ei[i] = 1.0
        x_phi = x0 + h*ei
        hessian[i,:] = (grad(x_phi) - g0)/h
    # hessian should be symmetric
    
    return hessian
    

def call_with_opts_from_dict(func):
    """
    call function func with its keywords replaced by values in dictionary.
    Avoids raising an error if dictionary contains keys that are not keywords
    in func.

    stolen from http://bytes.com/topic/python/answers/170128-expanding-dictionary-function-arguments
    """
    argnames = func.func_code.co_varnames[hasattr(func, 'im_func'):func.func_code.co_argcount]
    ndefaults = len(func.func_defaults or ())
    if ndefaults:
        optnames = argnames[-ndefaults:]
        argnames = argnames[:-ndefaults]
    else:
        optnames = []
    def _f(*args, **opts):
        try:
            actualargs = args #[args[argname] for argname in argnames]
            actualopts = {}
            for io,optname in enumerate(optnames):
                if optname in opts.keys(): 
                    actualopts[optname] = opts[optname]
                else:
                    actualopts[optname] = func.func_defaults[io]
#            print "Calling %s with arguments = %s and options = %s" % (func.func_name, actualargs, actualopts)
        except TypeError: raise TypeError, '%s(...) requires arg(s) %r'%(
            func.func_name, [argname for argname in argnames])
        return func(*actualargs, **actualopts)
    
    _f.func_name = func.func_name
    return _f

# Old parts of the code still import the option parser from this module.
# In order not to break this we have to import the option parser here.
from optparse import OptionParserFuncWrapper


if __name__ == "__main__":
    def bla(x,y,z=100,w=1000):
        print "x = %s, y = %s, z = %s, w = %s" % (x,y,z,w)
    opts = {'w': -1000, 'test':100}
    call_with_opts_from_dict(bla)(1,2,**opts)
