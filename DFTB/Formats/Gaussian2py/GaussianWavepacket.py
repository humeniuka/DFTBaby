import numpy as np
from numpy import linalg as la
from numpy import random
import time

def cw_to_AB(c,w):
    # c,w stand for center and width
    # w (q - c)^2 = w q^2 - 2 w * c * q + w * c^2
    A = -w
    B = -2*w*c
    return A,B

def AB_to_cw(A,B):
    w = A
    c = B/A
    return c,w

def gaussian_norm(A,B):
    """
    integral Int dQ g(Q) = Int dQ exp(-1/2 Q^T*A*Q + B^T*Q)
    """
    N = A.shape[0]
    Ainv = la.inv(A)
    nrm  = pow(2.0*np.pi,N/2.0)/np.sqrt(la.det(A)) 
    nrm *= np.exp(1.0/2.0 * np.dot(B, np.dot(Ainv, B)))
    return nrm
    

class Gaussian:
    """
    G(Q1,...,Q2^n) = exp( -1/2 Q^T*A*Q + B^T*Q )
    """
    def __init__(self, A, B):
        self.A = A
        self.B = B
    def __str__(self):
        txt = "Gaussian wavepacket: N*exp( -1/2 Q^T*A*Q + B^T*Q )\n"
        txt += " A:\n%s\n" % self.A
        txt += " B:\n%s\n" % self.B
        return txt
    def norm(self):
        return gaussian_norm(self.A, self.B)
    def evaluate(self, q):
        assert self.A.shape == (1,1)
        # exp(-1/2 q^T*A*q + B^T*q)
        psi = np.exp(-0.5 * np.dot(np.transpose(q), np.dot(self.A, q)) + np.dot(self.B, q))
        psi /= np.sqrt(gaussian_norm(2.0*self.A, 2.0*self.B))
        return psi
    def wigner_transform(self, hbar=1):
        print "Wigner transform"
        ReA = self.A.real
        ReAinv = la.inv(ReA)
        ImA = self.A.imag
        ReB = self.B.real
        ImB = self.B.imag
        Aw = np.bmat([[2*ReA + 2*np.dot(ImA, np.dot(ReAinv, ImA)),
                       2.0/hbar * np.dot(ImA, ReAinv)],
                      [2.0/hbar * np.dot(ReAinv, ImA), 
                       2.0/pow(hbar,2) * ReAinv]])
        Aw = np.asarray(Aw)
        Bw = 2 * np.hstack([ReB + np.dot(ImB, np.dot(ReAinv, ImA)), 
                        1.0/hbar * np.dot(ImB, ReAinv)])
        return Gaussian(Aw,Bw)
    def sample(self, nr_initial_conditions, pw=1.0):
        """
        Parameters:
        ===========
        nr_initial_conditions: number of samples that should be drawn
        pw: instead of sampling from the Wigner distribution W(q,p)
           the samples are drawn from the distribution [W(q,p)]^pw. 
           pw < 1.0 leads to a diffuser distribution of initial conditions
            

        Returns:
        ========
        qs, ps: Samples 
          qs[:,n], ps[:,n] are the initial positions and momenta of the n-th trajectory
        """
        # covariance matrix
#        cov = la.inv(self.A)
        cov = la.inv(pw*self.A)
        # mean
#        mu = np.dot(cov, self.B)
        mu = np.dot(cov, pw*self.B)
        # sample from multivariate normal distribution
        N = len(self.A)
        random.seed(int(time.time()*1000))
        samples = random.multivariate_normal(mu, cov, nr_initial_conditions).transpose()
        # dim = number of nuclear degrees of freedom
        dim = N/2
        qs = samples[:dim,:]
        ps = samples[dim:,:]
        # qs[i,:] -> initial positions on i-axis
        # ps[i,:] -> initial momenta along i-axis for all samples
        return (qs, ps)

if __name__ == "__main__":
    A = np.array([[1.0]])
    B = np.array([0.0])
    g = Gaussian(A,B)
    w = g.wigner_transform(hbar=1.0)
    q,p = w.sample(10)
    print "q = %s" % q
    print "p = %s" % p
