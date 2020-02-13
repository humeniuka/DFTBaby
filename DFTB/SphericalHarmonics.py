"""
implementation of spherical harmonics
"""
import numpy as np
from numpy import sqrt, pi
from scipy import special
from DFTB.MolecularIntegrals.SphericalCoords import cartesian2spherical


def Y(l,m,r):
	"""spherical harmonic Ylm in cartesian coordinates, r=[x,y,z]"""
	assert(l >= abs(m))
	(x,y,z) = (r[0], r[1], r[2])
        rn,th,phi = cartesian2spherical((x,y,z))
        return special.sph_harm(m,l,phi,th)

def Yold(l,m,r):
	"""first few spherical harmonic Ylm in cartesian coordinates, r=[x,y,z]"""
	assert(l >= abs(m))
	(x,y,z) = (r[0], r[1], r[2])
	rn = sqrt(x*x+y*y+z*z)
	if (l==0):
		if (m==0):
			return(1.0/2.0*sqrt(1.0/pi))		
	elif(l==1):
		if (m==-1):
			return(1.0/2.0*sqrt(3.0/(2.0*pi))*(x-1.0j*y)/rn)
		elif (m==0):
			return(1.0/2.0*sqrt(3.0/pi)*z/rn)
		elif (m==1):
			return(-1.0/2.0*sqrt(3.0/(2.0*pi))*(x+1.0j*y)/rn)
	elif(l==2):
		if (m==-2):
			return(1.0/4.0*sqrt(15.0/(2.0*pi))*pow( (x-1.0j*y)/rn, 2))
		elif(m==-1):
			return(1.0/2.0*sqrt(15.0/(2.0*pi))*(x-1.0j*y)*z/pow(rn,2))
		elif(m==0):
			return(1.0/4.0*sqrt(5.0/pi)*(2.0*z*z - x*x - y*y)/pow(rn,2))
		elif(m==1):
			return(-1.0/2.0*sqrt(15.0/(2.0*pi))*(x+1.0j*y)*z/pow(rn,2))
		elif(m==2):
			return(1.0/4.0*sqrt(15.0/(2.0*pi))*pow( (x+1.0j*y)/rn, 2))
	elif(l==3):
		if (m==-3):
			return (1.0/8.0*sqrt(35.0/pi) * pow( (x-1.0j*y)/rn, 3))
		elif (m==-2):
			return (1.0/4.0*sqrt(105.0/(2.0*pi)) * pow(x-1.0j*y,2) * z / pow(rn,3))
		elif (m==-1):
			return (1.0/8.0*sqrt(21.0/pi) * (x-1.0j*y)*(4.0*z*z - x*x - y*y)/pow(rn,3))
		elif (m==0):
			return (1.0/4.0*sqrt(7.0/pi)*z*(2.0*z*z-3.0*x*x-3.0*y*y)/pow(rn,3))
		elif (m==1):
			return (-1.0/8.0*sqrt(21.0/pi)*(x+1.0j*y)*(4.0*z*z-x*x-y*y)/pow(rn,3))
		elif (m==2):
			return (1.0/4.0*sqrt(105.0/(2.0*pi))*pow(x+1.0j*y,2)*z/pow(rn,3))
		elif (m==3):
			return (-1.0/8.0*sqrt(35.0/pi)*pow( (x+1.0j*y)/rn, 3))
        elif(l==4):
                if (m==-4):
                        return (3.0/16.0*sqrt(35.0/(2.0*pi))*pow( (x-1.0j*y)/rn, 4))
                elif (m==-3):
                        return (3.0/8.0*sqrt(35.0/pi)* (pow(x-1.0j*y,3)*z/pow(rn,4)))
                elif (m==-2):
                        return (3.0/8.0*sqrt(5.0/(2.0*pi))*pow(x-1.0j*y,2)*(7*z**2-rn**2)/pow(rn,4))
                elif (m==-1):
                        return (3.0/8.0*sqrt(5.0/pi)*(x-1.0j*y)*z*(7*z**2-3*rn**2)/pow(rn,4))
                elif (m==0):
                        return (3.0/16.0*sqrt(1.0/pi)*(35*z**4-30*z**2*rn**2+3*rn**4)/pow(rn,4))
                elif (m==1):
                        return (-3.0/8.0*sqrt(5.0/pi)*(x+1.0j*y)*z*(7*z**2-3*rn**2)/pow(rn,4))
                elif (m==2):
                        return (3.0/8.0*sqrt(5.0/(2.0*pi))*pow(x+1.0j*y,2)*(7*z**2-rn**2)/pow(rn,4))
                elif (m==3):
                        return (-3.0/8.0*sqrt(35.0/pi)*pow(x+1.0j*y,3)*z/pow(rn,4))
                elif (m==4):
                        return (3.0/16.0*sqrt(35.0/(2.0*pi))*pow((x+1.0j*y)/rn,4))                        
	else:
		raise Exception("Spherical harmonics for l>4 not implemented")

# real spherical harmonics
def Yreal(l,m,r):
	if m > 0:
		yreal = 1.0/sqrt(2.0) * (Y(l,m,r) + pow(-1.0,m)*Y(l,-m,r)) 
	elif m == 0:
		yreal = Y(l,m,r)
	else:
		# m < 0
		yreal = -1.0j/sqrt(2.0) * (Y(l,-m,r) - pow(-1.0,m)*Y(l,m,r))
	yreal *= pow(-1.0, m)
	return yreal.real


if __name__ == "__main__":
        x = np.random.rand(3)
        y = np.random.rand(3)
        z = np.random.rand(3)

        l = 3
        m = -2

        print Y(l,m,(x,y,z))
        print Yold(l,m,(x,y,z))
