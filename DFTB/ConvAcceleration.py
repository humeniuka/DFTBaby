"""
mixing scheme for accelerating convergence
in self consistent calculations by either interpolating
the density matrix or the KS hamiltonian
"""

############ DENSITY INTERPOLATION ###################
class DensityMixer:
    def next_approximation(self, P):
        pass
    def relative_change(self):
        pass

# DIIS ############################################################
"""
implementation of Pulay's DIIS algorithm to speed up iterative solution of
self charge consistent DFTB.

see http://vergil.chemistry.gatech.edu/notes/diis/node2.html
and Pulay (1980), "Convergence acceleration of iterative sequences. the case of SCF iteration"

"""
from numpy import zeros, ravel, reshape, dot, where, argmax
from numpy.linalg import solve, norm, LinAlgError
import scipy.linalg as la

class DIIS_80(DensityMixer):
    def __init__(self, m, relchange_start=0.5):
        """
        Parameters:
        ===========
        m: maximum number of trial vectors that are kept in memory
        relchange_start: DIIS is started only if the relative change
                         between steps is lower than this number
        """
        self.m = m
        self.residual_vectors = []
        self.trial_vectors = []
        self.start_flag = False
        self.relchange_start = relchange_start
    def nr_trial_vectors(self):
        return len(self.trial_vectors)
    def reset(self):
        self.residual_vectors = []
        self.trial_vectors = []        
        self.start_flag = False
    def add_trial_vector(self, p):
        """
        add the vector p to the set of iterative solutions
        """
        self.trial_vectors.append(p)
        if len(self.trial_vectors) >= 2:
            # guess error vectors from change
            # with respect to previous iteration
            p_i = self.trial_vectors[-2]
            p_iplus1 = self.trial_vectors[-1]
            dp_i = p_iplus1 - p_i
            self.residual_vectors.append(dp_i)
        # throw away the oldes trial vectors
        if len(self.trial_vectors) > self.m:
            self.trial_vectors.pop(0)
            self.residual_vectors.pop(0)
    def get_DIIS_approximate(self):
        """
        determine the best linear combination of the previous
        solution vectors.

        Returns:
        ========
        p_(i+1), dp_(i+1)
        where 
        p_(i+1) optimized solution vector, next residual vector dp_(i+1)  

        We can test for convergence by calculating norm(dp_(i+1))/norm(p_(i+1))
        """
        n = len(self.residual_vectors)
        assert n > 1
        A = zeros((n+1,n+1))
        for i in range(0, n):
            for j in range(0, n):
                A[i,j] = dot(self.residual_vectors[i], self.residual_vectors[j])
            A[i,n] = -1.0
            A[n,i] = -1.0
        A[n,n] = 0.0
        b = zeros(n+1)
        b[-1] = -1.0
        # solve A*c = b
        c = solve(A,b)
        assert abs(sum(c[:-1]) - 1.0) < 1.0e-10
        # p = sum_i c_i*p_i
        pDIIS = 0.0*self.trial_vectors[0]
        for i in range(0, n):
            pDIIS += c[i]*self.trial_vectors[i]
        # replace old trial vector by DIIS approximation
        self.trial_vectors[-1] = pDIIS
        self.residual_vectors[-1] = pDIIS - self.trial_vectors[-2]
        dp = self.residual_vectors[-1]
        return pDIIS, dp
    def next_approximation(self, p):
        """
        add the trial vector p and replace it with the DIIS approximate

        Parameters:
        ===========
        p: numpy array of arbitrary shape

        Returns:
        ========
        p_next: numpy array with same shape as p
        """
        shape = p.shape
        # add new trial vector
        pflat = ravel(p)
        self.add_trial_vector(pflat)
        if self.start_flag == True:
            try:
                pflat_next, dp = self.get_DIIS_approximate()
            except LinAlgError:
                print "DIIS failed: singular matrix!"
                pflat_next = 0.5*self.trial_vectors[-2] + 0.5*self.trial_vectors[-1]
                raise Exception("?")
            p_next = reshape(pflat_next, shape)
        else:
            if len(self.trial_vectors) > 3:
                relchange = norm(self.residual_vectors[-1])/norm(self.trial_vectors[-1])
                # start DIIS if density matrix does not change too much
                if relchange < self.relchange_start:
                    #print "start DIIS"
                    self.start_flag = True

            p_next = p
        return p_next
    def relative_change(self):
        """
        compute the relative change in the last iteration
          |p_(i+1) - p_i|/|p_i|
        """
        assert(self.m > 1)
        if len(self.residual_vectors) == 0:
            change = 1.0
        else:
            #change = norm(self.residual_vectors[-1])/norm(self.trial_vectors[-1])
            # average relative change over the last 2 iterations

            change = 0.0
            navg = min(len(self.residual_vectors), 2)
            for i in range(1, navg+1):
                change += norm(self.residual_vectors[-i])/norm(self.trial_vectors[-i])
            change /= float(navg)

            #
        return change


##### KS HAMILTONIAN INTERPOLATION ####################################

class FockMixer:
    def interpolate_hamiltonian(self, H):
        pass
    def set_error_vector(self, err):
        pass

class DIIS_82(FockMixer):
    """
    Pulay's "Improved SCF Convergence Acceleration" of 1982
    """
    def __init__(self, m, start_threshold=0.1):
        """
        Parameters:
        ===========
        m: maximum number of trial vectors that are kept in memory
        start_threshold: Initiate DIIS procedure if the largest error vector is smaller than this threshold
        """
        self.m = m
        self.start_threshold = start_threshold
        self.start_flag = False
        self.error_vectors = []
        self.error_norms = []
        self.trial_vectors = []
    def nr_trial_vectors(self):
        return len(self.trial_vectors)
    def reset(self):
        self.error_vectors = []
        self.error_norms = []
        self.trail_vectors = []
        self.start_flag = False
    def set_error_vector(self, err):
        err = ravel(err)
        self.error_vectors.append(err)
        self.error_norms.append(norm(err))
        assert len(self.error_vectors) == len(self.trial_vectors)
        # throw away the oldest trial vectors
        if len(self.trial_vectors) > self.m:
            self.trial_vectors.pop(0)
            self.error_vectors.pop(0)
            self.error_norms.pop(0)
    def get_DIIS_approximate(self):
        """
        determine the best linear combination of the previous
        solution vectors.

        Returns:
        ========
        v_(i+1) next optimized solution vector
        """
        n = len(self.error_vectors)
        assert n > 1
        A = zeros((n+1,n+1))
        for i in range(0, n):
            for j in range(0, n):
                A[i,j] = dot(self.error_vectors[i], self.error_vectors[j])
            A[i,n] = -1.0
            A[n,i] = -1.0
        A[n,n] = 0.0
        b = zeros(n+1)
        b[-1] = -1.0
        # solve A*c = b
        c = solve(A,b)
        print "best linear combination c = %s" % c
        assert abs(sum(c[:-1]) - 1.0) < 1.0e-5
        # v = sum_i c_i*v_i
        vDIIS = 0.0*self.trial_vectors[0]
        for i in range(0, n):
            vDIIS += c[i]*self.trial_vectors[i]
        return vDIIS
    def interpolate_hamiltonian(self, H):
        """
        Parameters:
        ===========
        H: current Kohn-Sham hamiltonian
        """
        H_next = H
        if self.start_flag == False and len(self.trial_vectors) > 2:
            # find the maximum relative error
            imax = argmax(self.error_norms)
            rel_err = self.error_norms[imax]/norm(self.trial_vectors[imax])
            print "threshold for DIIS: %s   relative error: %s" % (self.start_threshold, rel_err)
            if rel_err < self.start_threshold:
                print "start DIIS"
                self.start_flag = True

        if self.start_flag == True:
            try:
                vflat_next = self.get_DIIS_approximate()
            except LinAlgError:
                print "DIIS failed: singular matrix!"
                vflat_next = 0.5*self.trial_vectors[-2] + 0.5*self.trial_vectors[-1]
                raise Exception("?")
            H_next = reshape(vflat_next, H.shape)
        self.trial_vectors.append(ravel(H))
        return H_next
    def relative_change(self):
        """
        """
        assert(self.m > 1)
        if len(self.trial_vectors) < 2:
            change = 1.0
        else:
            change = norm(self.trial_vectors[-1] - self.trial_vectors[-2])/norm(self.trial_vectors[-1])
        return change
