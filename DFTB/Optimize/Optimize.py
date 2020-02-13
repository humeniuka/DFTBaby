# -*- coding: utf-8 -*-
"""
solve the following optimization problem:

  minimize f(x)      subject to  c_i(x) > 0   for  i=1,...,m

where f(x) is a scalar function, x is a real vector of size n, and c_i(x) are the 
m strict inequality constraints. The feasible space {x|C(x) > 0} is assumed to be convex.
The constraints are enforced by minimizing an auxiliary function f(x)+nu*B(x). 
B(x) is the log-barrier

  B(x) = - sum_i log(c_i(x))

and `nu` is a small adjustable number.

References
----------
 [1] J. Nocedal, S. Wright, 'Numerical Optimization', Springer, 2006
"""
import numpy as np
import numpy.linalg as la

from DFTB.Optimize import Constraints   # contains definition of log-barrier

from DFTB.utils import numerical_gradient, numerical_hessian_G

    
def line_search_backtracking(xk, fk, grad_fk, pk, func,
                             constraints=None,
                             a0=1.0, rho=0.3, c=0.0001, lmax=100):
    """
    perform a line search along the search direction pk using the Armijo backtracking algorithm

    see Algorithm 3.1 in Ref. [1]

    Parameters
    ----------
    xk: vector, current point
    fk: scalar, current value of objective function
    grad_fk: vector, current gradient 
    pk: vector, search direction
    func: callable, computes function value and gradient, 
          (fx, dfdx) = func(x)
    
    Optional
    --------
    constraints: callable, constraints(x) should return a vector `C` with the values 
              of the constraints at x and matrix `A` whose rows are the gradients 
              of the constraints w/r/t to x
    lmax: maximum number of tries before giving up

    The meaning of `a0`, `rho` and `c` can be glanced from the algorithm description in Ref.[1]

    Returns
    -------
    xkp1: next point
    """
    n = len(xk)
    a = a0
    # directional derivative
    df = np.dot(grad_fk, pk)
    # check that pk is a descent direction
    assert df <= 0.0, "pk=%s not a descent direction" % pk
    
    for l in range(0, lmax):
        x_interp = xk + a*pk
        # choose the step small enough so that no constraints are violated
        if constraints != None:
            C, A = constraints(x_interp)
            if np.any(C <= 0.0):
                # inequality constraints C > 0 were violated, a call to func() would
                # result in an error, reduce step size
                a *= rho
                continue
                        
        if func(x_interp) <= fk + c*a*df:
            break
        else:
            a *= rho
    else:
        raise RuntimeError("Linesearch failed! Could not find a step length that satisfies the sufficient decrease condition.")
    return x_interp

def line_search_wolfe(xk, fk, grad_fk, pk, func,
                      constraints=None, max_steplen=None, 
                      a0=1.0, amax=50.0, c1=0.0001, c2=0.9,
                      lmax=100, 
                      debug=0, debug_info=None):
    """
    find step size `a`` that satisfies the strong Wolfe conditions:  
                                                                      __
      1) sufficient decrease condition   f(xk + a pk) <= f(xk) + c1 a \/f(xk).pk
                                         __                        __
      2) curvature condition            |\/f(xk + a pk).pk| <= c2 |\/f(xk).pk|

    see Algorithm 3.5 and 3.6 in Ref. [1]

    Parameters
    ----------
    xk: vector, current point
    fk: scalar, current value of objective function
    grad_fk: vector, current gradient 
    pk: vector, search direction
    func: callable, computes function value and gradient, 
          (fx, dfdx) = func(x)
    
    Optional
    --------
    constraints: callable, constraints(x) should return a vector `C` with the values 
              of the constraints at x and matrix `A` whose rows are the gradients 
              of the constraints w/r/t to x
    max_steplen: callable, `max_steplen(x,v)` computes the largest step length `amax` such that
              x + amax*v  lies inside the feasible set or on its boundary.
    lmax: maximum number of tries before giving up

    For the meaning of `a0`, `amax`, `c1` and `c2` see the algorithm description in Ref.[1]

    Returns
    -------
    xkp1: next point
    """
    assert 0 < c1 < c2 < 1
    def s(a):
        """computes scalar function s: a -> f(xk + a*pk) and its derivative ds/da"""
        fx,dfdx = func(xk + a*pk)
        dsda = np.dot(dfdx, pk)
        return fx, dsda
    
    # 
    s0 = fk      # s(a=0.0)
    Ds0 = np.dot(grad_fk, pk)  # ds/da(a=0.0)
        
    def zoom(alo,ahi, slo,shi):
        """
        find a step length a that satisfies Wolfe's conditions by bisection inside in the interval [alo,ali] 
        """

        if debug > 0:
            print "alo=%e  ahi=%e   s(alo)=%e  s(ahi)=%e" % (alo, ahi, slo, shi)
            import matplotlib.pyplot as plt
            npts = 200
            a_arr = np.linspace(alo, ahi, npts)
            s_arr = [s(a)[0] for a in a_arr]
            #

            if debug_info != None:
                ensKin, ensCoul, ensXC = [], [], []
                for a in a_arr:
                    enKin, enCoul, enXC = debug_info(xk+a*pk)
                    ensKin.append( enKin )
                    ensCoul.append( enCoul )
                    ensXC.append( enXC )
                ensKin = np.array(ensKin)
                ensCoul = np.array(ensCoul)
                ensXC = np.array(ensXC)
                s_arr = np.array(s_arr)
                enKinMin = ensKin.min()
                enCoulMin = ensCoul.min()
                enXCMin = ensXC.min()
                plt.plot(a_arr, ensKin-enKinMin, label="kinetic")
                plt.plot(a_arr, ensCoul-enCoulMin, label="Coulomb")
                plt.plot(a_arr, ensXC-enXCMin, label="exchange-correlation")
                plt.plot(a_arr, s_arr-s_arr.min(), label="$s(\\alpha)$")

            else:
                plt.plot(a_arr, s_arr, label="$s(\\alpha)$")

            #plt.plot(a_arr, s_arr, label="$s(\\alpha)$")
            plt.legend()
            plt.show()

        for j in range(0, lmax):
            # evaluate s and s' at the midpoint of the search interval
            aj = 0.5*(alo+ahi)
            sj,Dsj = s(aj)
            
            if debug > 0:
                print "aj=%e  alo=%e  ahi=%e" % (aj, alo, ahi)
                plt.xlabel("step length $\\alpha$")
                plt.plot(a_arr, s_arr, label="$s(\\alpha)$")
                plt.plot([ahi],[shi], "^", color="red")
                plt.plot([alo],[slo], "v", color="red")
                plt.plot([aj], [sj],  "x", color="red")
                # tangent at s(aj)
                plt.plot([aj,ahi],[sj,sj+Dsj*(ahi-aj)], color="green", label="$s'(\\alpha_j)$")
                plt.plot([alo,aj],[sj+Dsj*(alo-aj),sj], color="green")

                plt.show()

            if (sj > s0 + c1*aj*Ds0) or (sj >= slo):
                # sufficient decrease condition is not fulfilled
                ahi = aj
                shi = sj
            else:
                if abs(Dsj) <= c2*abs(Ds0):
                    # curvature condition met, we are done
                    aWolfe = aj
                    break
                if Dsj*(ahi-alo) >= 0.0:
                    ahi = alo
                    shi = slo
                alo = aj
                slo = sj
        else:
            msg = "``zoom`` could not find a point satisfying Wolfe's condition in the interval [%e,%e] in %d iterations!" % (alo,ahi, j)
            aWolfe = 0.5*(alo+ahi)
            print "WARNING: %s" % msg
            #raise RuntimeError(msg)
            
        return aWolfe
            
    def feasible(a):
        """checks whether the point xk+a*pk is inside the feasible region"""
        if constraints != None:
            C,A = constraints(xk + a*pk)
            if np.any(C <= 0.0):
                # inequality constraints C > 0 are violated
                return False
        return True
    
    # Find the largest feasible step length.
    if max_steplen != None:
        _amax = max_steplen(xk, pk)
        if _amax == np.inf:
            # If there is no limit on the maximum step length, we leave the default value `amax` untouched
            pass
        else:
            # override default value with maximum step length
            amax = _amax
            ####
            """
            if debug > 0:
                feasible_arr = []
                amax_arr = np.linspace(0.0, 10*amax, 2000)
                for _amax in amax_arr:
                    if feasible(_amax):
                        feasible_arr.append(1)
                    else:
                        feasible_arr.append(0)
                if not np.all(feasible_arr):
                    import matplotlib.pyplot as plt
                    plt.xlabel("amax")
                    plt.plot(amax_arr, feasible_arr)
                    plt.show()
            """
            ####
#            # any step length < `amax` should lie in the feasible region,
#            # while any step length > `amax` should lie outside
#            assert feasible(0.999*amax) == True
#            assert feasible(1.001*amax) == False
    # Since the feasible region is convex, all step lengths `a` < `amax`
    # will be feasible as well, if amax is feasible.
    else:
        # If no function was provided for selecting the maximum step length,
        # we have to ensure that `amax` at least does not lie outside the feasible set.
        while feasible(amax) == False:
            # decrease `amax` until it lies in the feasible region
            amax *= 0.99
    # The initial guess for the step length should satisfy `a0` < `amax`.
    if a0 >= amax:
        a0 = 0.5*amax

    if debug > 1:
        import matplotlib.pyplot as plt
        npts = 200
        a_arr = np.linspace(0.0, amax, npts)
        s_arr = [s(a)[0] for a in a_arr]
        plt.plot(a_arr, s_arr, label="$s(\\alpha)$")
        plt.legend()
        plt.show()
        
    #print "amax=%s   a0=%s" % (amax, a0)
    # Algorithm 3.5, brackets the interval 
    aim1 = 0.0
    sim1 = s0
    ai = a0
    for i in range(1, lmax):
        si,Dsi = s(ai)
        if (si > s0 + c1*ai*Ds0) or ((si >= sim1) and i > 1):
            # sufficient decrease condition is not fulfilled => a minimum has to lie
            # in between around, which the Wolfe conditions hold.
            aWolfe = zoom(aim1,ai, sim1,si)
            break
        if abs(Dsi) <= c2*abs(Ds0):
            # curvature condition is fulfilled
            aWolfe = ai
            break
        if Dsi >= 0.0:
            # derivative s'(a) changed sign => a minimum has to lie in between
            aWolfe = zoom(ai,aim1, si,sim1)
            break
        aim1 = ai
        sim1 = si
        # choose a new a_(i+1) from the interval (ai,amax)
        ai = 2*aim1
        if ai >= amax:
            print "WARNING: end of search interval reached!"
            # end of interval reached when multiplying ai by 2, approach amax from below
            ai = 0.5*(aim1 + amax)
        #print "a_(i-1) = %e  a_i = %e" % (aim1, ai)
    else:
        raise RuntimeError("Linesearch failed! Could not find a step length that satisfies Wolfe's conditions.")
    x_Wolfe = xk + aWolfe*pk
    return x_Wolfe
            
def bfgs_update(invHk, sk, yk, k):
    """
    update the inverse Hessian invH_(k+1) based on Algorithm 6.1 in Ref.[1]

    Parameters:
    -----------
    invHk: inverse Hessian approximation __            __
    yk: gradient difference vector, yk=  \/ f_(k+1)  - \/ f_k
    sk: step vector,  x_(k+1) - x_k
    k: integer counting the number of iterations, 
        if k=0, 

    Returns
    -------
    invH_(k+1): next inverse Hessian approximation
    """
    n = len(sk)
    Id = np.eye(n)
    assert k >= 1
    if k == 1:
        invHkp1 = np.dot(yk,sk)/np.dot(yk,yk) * Id
    else:
        rk = 1.0/np.dot(yk,sk)
        U = Id - rk*np.outer(sk,yk)
        V = Id - rk*np.outer(yk,sk)
        W = rk*np.outer(sk,sk)
        
        invHkp1 = np.dot(U, np.dot(invHk, V)) + W
    return invHkp1
    
class OptimizationResult:
    def __init__(self, x, fun, grad, nit):
        self.x = x
        self.fun = fun
        self.grad = grad
        self.nit = nit


def minimize(objfunc, x0,
             method="BFGS",
             line_search_method="Wolfe", 
             constraints=None, max_steplen=None, bounds=None,
             callback=None, maxiter=100000,
             gtol=1.0e-6, ftol=1.0e-8,
             debug=0, debug_info=None):
    """
    minimize a scalar function ``objfunc``(x) possibly subject to constraints.

    The minimization is converged if
      * |df/dx| < gtol and 
      * |f(k+1)-f(k)| < ftol

    Parameters
    ----------
    objfunc: callable, objective function that should be minimized, it should return
             the function value and gradient:

               fx, dfdx = objfunc(x) 

    x0: initial point, where the optimization starts

    Optional
    --------
    method: choose how the search direction should be determined,
       'Newton'           - the Hessian H is calculated by numerically differentiating the gradient,
                            the search direction is then obtained by solving  H.p = -df/dx
       'Steepest Descent' - the search direction is antiparallel to the gradient, p = -df/dx
       'BFGS'             - an approximation to the inverse Hessian invH is updated after each step
                            so as to track the curvature along the path, then p = -invH.df/dx
    line_search_method: choose how the step length in the search direction should be determined
       'Armijo'           - starting with a=1, the algorithm backtracks until sufficient decrease is obtained
       'Wolfe'            - attempts to find a step length that satisfies Wolfe's conditions
       'largest'          - takes the largest step compatible with max_steplen() along the search direction
    constraints: callable, constraints(x) should return a vector C with the values 
              of the constraints at x and matrix A whose rows are the gradients 
              of the constraints w/r/t to x
    max_steplen: callable, `max_steplen(x,v)` computes the largest step length `amax` such that
              x + amax*v  lies inside the feasible set or on its boundary.
    bounds: list of tuples (lower,upper) with bounds for each component of the vector x,
    callback: callable, at the end of each iteration this function is called with the current vector x
              as argument
    maxiter: maximum number of iterations
    gtol: tolerance for norm of gradient
    ftol: tolerance for change of function value

    Returns
    -------
    res: instance of OptimizationResult

    Notes
    -----
    The 'BFGS' algorithm should be combined with a 'Wolfe' line search. The curvature condition is important
    in a quasi-Newton method, because it ensures that the approximation to the Hessian remains positive definite
    after each update.

    See Also
    --------
    test_newton
    """    
    assert method in ["Newton", "Steepest Descent", "BFGS"]
    assert line_search_method in ["Armijo", "Wolfe", "largest"]
    n = len(x0)
    def barrier(x):
        nu = 1.0e-5  #0.001 #0.00000001 #0.001
        if constraints != None:
            C,A = constraints(x)
            # add log-barrier
            Bx,dBdx = Constraints.log_barrier(C,A)
        else:
            Bx = 0.0
            dBdx = np.zeros(n)
        Bx *= nu
        dBdx *= nu
        return Bx, dBdx
    def func(x):
        fx,dfdx = objfunc(x)
        Bx,dBdx = barrier(x)
        return fx+Bx
    def grad(x):
        fx,dfdx = objfunc(x)
        Bx,dBdx = barrier(x)
        return dfdx+dBdx
    def func_grad(x):
        fx,dfdx = objfunc(x)
        Bx,dBdx = barrier(x)
        return fx+Bx, dfdx+dBdx
    def hess(x):
        H = numerical_hessian_G(grad, x)
        return H

    if debug > 0:
        print constraints
        C = Contours(x0, func, constraints=constraints, bounds=bounds, xind=0, yind=1)

    xk = x0
    fk, grad_fk = func_grad(xk)
    converged = False
    # smallest representable positive number such that 1.0+eps != 1.0.
    epsilon = np.finfo(float).eps
    for k in range(0, maxiter):
        # determine new search direction
        if method == "Newton":
            # compute exact hessian numerically
            Hk = hess(xk)
            # make Hk sufficiently positive definite
            Bk, tau_k = modified_cholesky(Hk)
            # search direction
            pk = la.solve(Bk, -grad_fk)
        elif method == "Steepest Descent":
            pk = -grad_fk    # steepest descent
        elif method == "BFGS":
            if k == 0:
                invHk = np.eye(n)
            else:
                #assert np.dot(yk,sk) > 0.0
                if np.dot(yk,sk) <= 0.0:
                    print "WARNING: positive definiteness of Hessian approximation lost in BFGS update, since yk.sk <= 0!"
                #
                invHk = bfgs_update(invHk, sk, yk, k)
            pk = np.dot(invHk,-grad_fk)
        # 
        if debug > 0:
            C.plot_contours(xk, fk, grad_fk, pk, Bk)
        # determine next point by a line search
        if line_search_method == "Armijo":
            x_kp1 = line_search_backtracking(xk, fk, grad_fk, pk, func, constraints=constraints)
        elif line_search_method == "Wolfe":
            # Quasi-Newton methods may fail to work well with backtracking line search, since
            # the curvature condition might not be fulfilled. In this case tiny steps are taken
            # along the descent direction, although a single larger step could be taken.
            ### DEBUG
            if k >= 3:
                _debug = 1
            else:
                _debug = 0
            _debug = 0
            ###
            x_kp1 = line_search_wolfe(xk, fk, grad_fk, pk, func_grad,
                                      constraints=constraints, max_steplen=max_steplen,
                                      debug=_debug, debug_info=debug_info)
        elif line_search_method == "largest":
            # simply take a large step along the search direction
            if max_steplen != None:
                amax = max_steplen(xk,pk)
            else:
                amax = 1.0
            x_kp1 = xk + amax*pk
            
        f_kp1, grad_f_kp1 = func_grad(x_kp1)
        # compute change of function value from step k to the next and norm of the gradient
        f_change = abs(f_kp1 - fk)
        gnorm = la.norm(grad_f_kp1)
        if f_change < ftol and gnorm < gtol:
            converged = True
        if f_change < epsilon:
            # f(k+1) and f(k) cannot be distinguished properly because of finite numerical precision
            print "WARNING: |f(k+1) - f(k)| < epsilon  (numerical precision) !"
            converged = True
        # step vector
        sk = x_kp1 - xk
        # gradient difference vector
        yk = grad_f_kp1 - grad_fk
        # new variables for step k become old ones for step k+1
        xk = x_kp1
        fk = f_kp1
        grad_fk = grad_f_kp1
        if callback != None:
            callback(xk)
        print "k=%10.1d  f(x) = %15.10f  |x(k+1)-x(k)| = %e  |f(k+1)-f(k)| = %e  |df/dx| = %e" % (k, fk, la.norm(sk), f_change, gnorm)
        if converged == True:
            #if method == "Newton":
            #    # tau is the smallest positive number, such that H+tau*Id is positive definite
            #    assert abs(tau_k) == 0.0, "At the minimum the Hessian should be positive definite. (i.e. tau_k = 0), but got tau_k = %s!" % tau_k
            break
    else:
        raise RuntimeError("No convergence in Newton's method after %d iterations!" % (k+1))
    return OptimizationResult(xk, fk, grad_fk, k)
    
def modified_cholesky(A, beta=0.01):
    """
    make the matrix A positive definite by adding a small positive multiple of the identity:
      A -> A + tau*Id

    see Algorithm 3.3 in Ref. [1]
    """
    n,n = A.shape
    
    min_diag = np.diag(A).min()
    if min_diag > 0.0:
        tau_0 = 0.0
    else:
        tau_0 = -min_diag + beta
    tau_k = tau_0
    #print "tau_0 = %e" % tau_0
    Id = np.eye(n)
    while True:
        #print "  tau_k = %e" % tau_k 
        try:
            L = la.cholesky(A + tau_k * Id)
        except la.LinAlgError as e:
            #print "tau_k = %s   %s" % (tau_k, str(e))
            tau_k = max(2*tau_k, beta)
            continue
        break
    #print "tau_k = %e" % tau_k
    Apos = A + tau_k*Id
    return Apos, tau_k

############ for TESTING ########

import copy
#import matplotlib.pyplot as plt

class Contours:
    """
    compare the contours of func(x) with the contours of the quadratic model
       mk(x) = func(xk) + g^T.x + 1/2 x^T.B.x
    """
    def __init__(self, point_0, func, constraints=None, bounds=None, xind=0, yind=1):
        self.xind = xind
        self.yind = yind
        self.points = [point_0]   # list of the points visited during the optimization
        self.X, self.Y, self.Z = self.eval_function(point_0, func,
                                                    constraints=constraints, bounds=bounds,
                                                    xind=self.xind, yind=self.yind)
    def plot_contours(self, point_k, f, g, p, B):
        print "current point = %s" % point_k
        print "current function value = %s" % f
        print "current gradient = %s" % g
        print "current search direction = %s" % p
        self.points.append( point_k )
        model = self.quadratic_model(point_k, f, g, B, np.array((self.X, self.Y)))
        CS_objfunc = plt.contour(self.X, self.Y, self.Z, 50, colors="black")
        CS_model = plt.contour(self.X, self.Y, model, 40, colors="red")
        plt.clabel(CS_objfunc, inline=1)
        plt.clabel(CS_model, inline=1)
        #
        x_pos = []
        y_pos = []
        for point_k in self.points:
            x_pos.append(point_k[self.xind])
            y_pos.append(point_k[self.yind])
                         
        plt.plot(x_pos, y_pos, color="black")
        plt.plot(x_pos, y_pos, "o", color="red")
        # show search direction
        ax = plt.axes()
        ax.arrow(point_k[self.xind], point_k[self.yind], p[self.xind], p[self.yind])
        
        plt.show()
        
    def eval_function(self, point_k, func, constraints=None, bounds=None, xind=0, yind=1):
        d = 0.9
        x0 = point_k[xind]
        y0 = point_k[yind]
        if bounds[xind][0] != None:
            # only width parameter has lower bound  a >= 0
            xmin, xmax = (1-d)*x0, (1+d)*x0
        else:
            xmin, xmax = -(1+d)*x0, (1+d)*x0
        if bounds[yind][0] != None:
            # only width parameter has lower bound  a >= 0
            ymin, ymax = (1-d)*y0, (1+d)*y0
        else:
            ymin, ymax = -(1+d)*y0, (1+d)*y0
            
        N = 150
        X, Y = np.mgrid[xmin:xmax:N*1j, ymin:ymax:N*1j]
        shape = X.shape
        point_i = copy.deepcopy(point_k)
        zvec = []
        for xi,yi in zip(X.flatten(), Y.flatten()):
            point_i[xind] = xi
            point_i[yind] = yi
            if constraints != None:
                C, A = constraints( point_i )
                if np.any(C <= 0.0):
                    # inequality constraints C > 0 were violated, a call to func() would
                    # result in an error
                    fx = np.nan
                else:
                    fx = func( point_i )
            else:
                fx = func( point_i )
            zvec.append(fx)
        zvec = np.array(zvec)
        Z = zvec.reshape(shape)
        
        return X,Y,Z
    def quadratic_model(self, xk, f, g, B,  x):
        """
        evaluate the quadratic approximation m(x) on a grid x
        """
        m = f
        indeces = [self.xind, self.yind]
        for i in range(0, 2):
            m += g[i]*(x[i,...]-xk[i])
            for j in range(0, 2):
                m += 0.5*(x[i,...]-xk[i])*B[indeces[i],indeces[j]]*(x[j,...]-xk[j])
        return m

def check_gradient(f, x0):
    """compare numerical and analytical gradient"""
    from scipy.optimize import check_grad
    def func(x):
        return f(x)[0]
    def grad(x):
        return f(x)[1]
    err_rel = check_grad(func, grad, x0)/func(x0)
    print "grad (NUMERICAL)  = %s" % numerical_gradient(func, x0)
    print "grad (ANALYTICAL) = %s" % grad(x0)
    assert err_rel < 1.0e-5, "relative error of gradient = %e !" % err_rel
    
    return err_rel
    
############ (constrained) Rosenbrock's function ###########################
    
def rosenbrock(x, a=1.0, b=100.0):
    f = (a-x[0])**2 + b*(x[1]-x[0]**2)**2
    dfdx = 0*x
    dfdx[0] = -2*(a-x[0]) - 4*b*x[0]*(x[1]-x[0]**2)
    dfdx[1] = 2*b*(x[1]-x[0]**2)
    
    return f, dfdx

def rosenbrock_constraints(x, a=1.0, b=100.0):
    """
    constraint   
       x1^2 + x2^2 < 1
    and its gradients
    """
    c1 = 1.0 - x[0]**2 - x[1]**2
    C = np.array([c1])
    A = np.array([[-2*x[0], -2*x[1]]])
    return C, A
    
def test_newton():
    x0 = np.array([0.0, 0.0])
    ### check gradient
    check_gradient(rosenbrock, x0)
    ### check minimization routine

    # minimize Rosenbrock's function without constraints
    res = minimize(rosenbrock, x0)
    print res.x, res.fun
    # minimize Rosenbrock's function subject to  c1 = 1 - x1^2 - x2^2 > 0
    # ... with Newton (using exact Hessian)
    print "Newton with Armijo line search"
    res = minimize(rosenbrock, x0, method="Newton", line_search_method="Armijo", constraints=rosenbrock_constraints)
    print res.x, res.fun
    print "Newton with Wolfe line search"
    res = minimize(rosenbrock, x0, method="Newton", line_search_method="Wolfe", constraints=rosenbrock_constraints)
    print res.x, res.fun
    # ... with BFGS approximation to inverse Hessian
    print "BFGS with Armijo line search"
    res = minimize(rosenbrock, x0, method="BFGS", line_search_method="Armijo", constraints=rosenbrock_constraints)
    print res.x, res.fun
    print "BFGS with Wolfe line search"
    res = minimize(rosenbrock, x0, method="BFGS", line_search_method="Wolfe", constraints=rosenbrock_constraints)
    print res.x, res.fun
    # ... with Steepest Descent
    print "Steepest Descent"
    res = minimize(rosenbrock, x0, method="Steepest Descent", constraints=rosenbrock_constraints, gtol=1.0e-7)
    print res.x, res.fun
    
if __name__ == "__main__":
    test_newton()

    
