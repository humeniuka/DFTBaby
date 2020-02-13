import numpy as np
from scipy import interpolate
import _pylibxc

def list_functionals():
    for func_id in range(1, 245):
        try:
            func = _pylibxc.libXCFunctional(func_id)
        except ValueError:
            pass
        print "================================="
        print func.description()

def interpolated_derivative(x, f):
    # interpolate f on a smooth grid and
    # compute gradient at positions x
    srt = np.argsort(x)
    xg,fg = x[srt], f[srt]
    tck = interpolate.splrep(xg, fg)
    dfdx = interpolate.splev(x, tck, der=1)
    return dfdx

class libXCFunctional:
    def __init__(self, exchange_name, correlation_name):
        self.func_x = _pylibxc.libXCFunctional(exchange_name)
        self.func_c = _pylibxc.libXCFunctional(correlation_name)
    def setGrid(self, x):
        self.x = x
    def getSigma(self, rho):
        if hasattr(self, "x"):
            # gradient in spherical coordinates points along e_r
            dndr = interpolated_derivative(self.x, rho)
            sigma = dndr*dndr
        else:
            sigma = 0.0*rho
        return sigma
    def exc(self, rho):
        sigma = self.getSigma(rho)
        ex = self.func_x.exc(rho, sigma)
        ec = self.func_c.exc(rho, sigma)
        return ex+ec
    def vxc(self, rho):
        sigma = self.getSigma(rho)
        vx = self.func_x.vxc(rho, sigma)
        vc = self.func_c.vxc(rho, sigma)
        return vx+vc
    def fxc(self, rho):
        sigma = self.getSigma(rho)
        fx = self.func_x.fxc(rho, sigma)
        fc = self.func_c.fxc(rho, sigma)
        return fx+fc

    def __str__(self):
        txt  = "Exchange\n"
        txt += "========\n"
        txt += self.func_x.description() + "\n"
        txt += "Correlation\n"
        txt += "===========\n"
        txt += self.func_c.description() + "\n"

if __name__ == "__main__":
    list_functionals()
