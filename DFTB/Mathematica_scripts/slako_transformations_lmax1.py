# This file has been generated automatically
from numpy import sin,cos,sqrt
phi = [ \
	lambda th1,th2: 0.5, \
	lambda th1,th2: 0, \
	lambda th1,th2: (cos(th2)*sqrt(3))/2., \
	lambda th1,th2: (cos(th1)*sqrt(3))/2., \
	lambda th1,th2: (3*sin(th1)*sin(th2))/4., \
	lambda th1,th2: (3*cos(th1)*cos(th2))/2., \
]

angular_phi = {\
	(0,0,0,0) : phi[0],	\
	(0,0,1,0) : phi[2],	\
	(1,0,0,0) : phi[3],	\
	(1,-1,1,-1) : phi[4],	\
	(1,0,1,0) : phi[5],	\
	(1,1,1,1) : phi[4],	\
}

tau2index = {\
	(0,0,0,0) : 0,\
	(0,0,1,0) : 2,\
	(1,0,0,0) : 3,\
	(1,-1,1,-1) : 4,\
	(1,0,1,0) : 5,\
	(1,1,1,1) : 4,\
}

index2tau = {\
	0 : (0,0,0,0),\
	2 : (0,0,1,0),\
	3 : (1,0,0,0),\
	4 : (1,-1,1,-1),\
	5 : (1,0,1,0),\
	4 : (1,1,1,1),\
}

# transformation rules for matrix elements
# x,y,z are directional cosines, r is the distance between the two centers
slako_transformations = {\
	(0,0,0,0) : lambda r,x,y,z,SorH: SorH(0,r), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: y*SorH(2,r), \
	(0,0,1,0) : lambda r,x,y,z,SorH: z*SorH(2,r), \
	(0,0,1,1) : lambda r,x,y,z,SorH: x*SorH(2,r), \
	(1,-1,0,0) : lambda r,x,y,z,SorH: y*SorH(3,r), \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (x**2 + z**2)*SorH(4,r) + y**2*SorH(5,r), \
	(1,-1,1,0) : lambda r,x,y,z,SorH: y*z*(-SorH(4,r) + SorH(5,r)), \
	(1,-1,1,1) : lambda r,x,y,z,SorH: x*y*(-SorH(4,r) + SorH(5,r)), \
	(1,0,0,0) : lambda r,x,y,z,SorH: z*SorH(3,r), \
	(1,0,1,-1) : lambda r,x,y,z,SorH: y*z*(-SorH(4,r) + SorH(5,r)), \
	(1,0,1,0) : lambda r,x,y,z,SorH: (x**2 + y**2)*SorH(4,r) + z**2*SorH(5,r), \
	(1,0,1,1) : lambda r,x,y,z,SorH: x*z*(-SorH(4,r) + SorH(5,r)), \
	(1,1,0,0) : lambda r,x,y,z,SorH: x*SorH(3,r), \
	(1,1,1,-1) : lambda r,x,y,z,SorH: x*y*(-SorH(4,r) + SorH(5,r)), \
	(1,1,1,0) : lambda r,x,y,z,SorH: x*z*(-SorH(4,r) + SorH(5,r)), \
	(1,1,1,1) : lambda r,x,y,z,SorH: (y**2 + z**2)*SorH(4,r) + x**2*SorH(5,r), \
}
# transformation rules for gradients of matrix elements
slako_transformations_gradient = [\
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: x*SorH(0,r,1), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: x*y*(-(SorH(2,r)/r) + SorH(2,r,1)), \
	(0,0,1,0) : lambda r,x,y,z,SorH: x*z*(-(SorH(2,r)/r) + SorH(2,r,1)), \
	(0,0,1,1) : lambda r,x,y,z,SorH: -(((-1 + x**2)*SorH(2,r))/r) + x**2*SorH(2,r,1), \
	(1,-1,0,0) : lambda r,x,y,z,SorH: x*y*(-(SorH(3,r)/r) + SorH(3,r,1)), \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (x*(-2*(-1 + x**2 + z**2)*SorH(4,r) + r*(x**2 + z**2)*SorH(4,r,1) + y**2*(-2*SorH(5,r) + r*SorH(5,r,1))))/r, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: (x*y*z*(2*SorH(4,r) - 2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1))))/r, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: (y*((-1 + 2*x**2)*SorH(4,r) + SorH(5,r) + x**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,0,0,0) : lambda r,x,y,z,SorH: x*z*(-(SorH(3,r)/r) + SorH(3,r,1)), \
	(1,0,1,-1) : lambda r,x,y,z,SorH: (x*y*z*(2*SorH(4,r) - 2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1))))/r, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (x*(-2*(-1 + x**2 + y**2)*SorH(4,r) + r*(x**2 + y**2)*SorH(4,r,1) + z**2*(-2*SorH(5,r) + r*SorH(5,r,1))))/r, \
	(1,0,1,1) : lambda r,x,y,z,SorH: (z*((-1 + 2*x**2)*SorH(4,r) + SorH(5,r) + x**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,1,0,0) : lambda r,x,y,z,SorH: -(((-1 + x**2)*SorH(3,r))/r) + x**2*SorH(3,r,1), \
	(1,1,1,-1) : lambda r,x,y,z,SorH: (y*((-1 + 2*x**2)*SorH(4,r) + SorH(5,r) + x**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,1,1,0) : lambda r,x,y,z,SorH: (z*((-1 + 2*x**2)*SorH(4,r) + SorH(5,r) + x**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (x*(-2*(-1 + x**2)*SorH(5,r) - (y**2 + z**2)*(2*SorH(4,r) - r*SorH(4,r,1)) + r*x**2*SorH(5,r,1)))/r, \
},
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: y*SorH(0,r,1), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: -(((-1 + y**2)*SorH(2,r))/r) + y**2*SorH(2,r,1), \
	(0,0,1,0) : lambda r,x,y,z,SorH: y*z*(-(SorH(2,r)/r) + SorH(2,r,1)), \
	(0,0,1,1) : lambda r,x,y,z,SorH: x*y*(-(SorH(2,r)/r) + SorH(2,r,1)), \
	(1,-1,0,0) : lambda r,x,y,z,SorH: -(((-1 + y**2)*SorH(3,r))/r) + y**2*SorH(3,r,1), \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (y*(-2*(-1 + y**2)*SorH(5,r) - (x**2 + z**2)*(2*SorH(4,r) - r*SorH(4,r,1)) + r*y**2*SorH(5,r,1)))/r, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: (z*((-1 + 2*y**2)*SorH(4,r) + SorH(5,r) + y**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: (x*((-1 + 2*y**2)*SorH(4,r) + SorH(5,r) + y**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,0,0,0) : lambda r,x,y,z,SorH: y*z*(-(SorH(3,r)/r) + SorH(3,r,1)), \
	(1,0,1,-1) : lambda r,x,y,z,SorH: (z*((-1 + 2*y**2)*SorH(4,r) + SorH(5,r) + y**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (y*(-2*(-1 + x**2 + y**2)*SorH(4,r) + r*(x**2 + y**2)*SorH(4,r,1) + z**2*(-2*SorH(5,r) + r*SorH(5,r,1))))/r, \
	(1,0,1,1) : lambda r,x,y,z,SorH: (x*y*z*(2*SorH(4,r) - 2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1))))/r, \
	(1,1,0,0) : lambda r,x,y,z,SorH: x*y*(-(SorH(3,r)/r) + SorH(3,r,1)), \
	(1,1,1,-1) : lambda r,x,y,z,SorH: (x*((-1 + 2*y**2)*SorH(4,r) + SorH(5,r) + y**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,1,1,0) : lambda r,x,y,z,SorH: (x*y*z*(2*SorH(4,r) - 2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1))))/r, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (y*(-2*(-1 + y**2 + z**2)*SorH(4,r) + r*(y**2 + z**2)*SorH(4,r,1) + x**2*(-2*SorH(5,r) + r*SorH(5,r,1))))/r, \
},
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: z*SorH(0,r,1), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: y*z*(-(SorH(2,r)/r) + SorH(2,r,1)), \
	(0,0,1,0) : lambda r,x,y,z,SorH: -(((-1 + z**2)*SorH(2,r))/r) + z**2*SorH(2,r,1), \
	(0,0,1,1) : lambda r,x,y,z,SorH: x*z*(-(SorH(2,r)/r) + SorH(2,r,1)), \
	(1,-1,0,0) : lambda r,x,y,z,SorH: y*z*(-(SorH(3,r)/r) + SorH(3,r,1)), \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (z*(-2*(-1 + x**2 + z**2)*SorH(4,r) + r*(x**2 + z**2)*SorH(4,r,1) + y**2*(-2*SorH(5,r) + r*SorH(5,r,1))))/r, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: (y*((-1 + 2*z**2)*SorH(4,r) + SorH(5,r) + z**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: (x*y*z*(2*SorH(4,r) - 2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1))))/r, \
	(1,0,0,0) : lambda r,x,y,z,SorH: -(((-1 + z**2)*SorH(3,r))/r) + z**2*SorH(3,r,1), \
	(1,0,1,-1) : lambda r,x,y,z,SorH: (y*((-1 + 2*z**2)*SorH(4,r) + SorH(5,r) + z**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (z*(-2*(-1 + z**2)*SorH(5,r) - (x**2 + y**2)*(2*SorH(4,r) - r*SorH(4,r,1)) + r*z**2*SorH(5,r,1)))/r, \
	(1,0,1,1) : lambda r,x,y,z,SorH: (x*((-1 + 2*z**2)*SorH(4,r) + SorH(5,r) + z**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,1,0,0) : lambda r,x,y,z,SorH: x*z*(-(SorH(3,r)/r) + SorH(3,r,1)), \
	(1,1,1,-1) : lambda r,x,y,z,SorH: (x*y*z*(2*SorH(4,r) - 2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1))))/r, \
	(1,1,1,0) : lambda r,x,y,z,SorH: (x*((-1 + 2*z**2)*SorH(4,r) + SorH(5,r) + z**2*(-2*SorH(5,r) + r*(-SorH(4,r,1) + SorH(5,r,1)))))/r, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (z*(-2*(-1 + y**2 + z**2)*SorH(4,r) + r*(y**2 + z**2)*SorH(4,r,1) + x**2*(-2*SorH(5,r) + r*SorH(5,r,1))))/r, \
},
]
slako_transformations_hessian = [\
[\
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: -(((-1 + x**2)*SorH(0,r,1))/r) + x**2*SorH(0,r,2), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*y*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x**2*y*SorH(2,r,2), \
	(0,0,1,0) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*z*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x**2*z*SorH(2,r,2), \
	(0,0,1,1) : lambda r,x,y,z,SorH: (3*x*(-1 + x**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x**3*SorH(2,r,2), \
	(1,-1,0,0) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*y*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x**2*y*SorH(3,r,2), \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (2*(-1 + 4*x**2)*(-1 + x**2 + z**2)*SorH(4,r) + 2*(-1 + 4*x**2)*y**2*SorH(5,r) + r*((z**2 - 5*x**2*(-1 + x**2 + z**2))*SorH(4,r,1) + (1 - 5*x**2)*y**2*SorH(5,r,1) + r*x**2*((x**2 + z**2)*SorH(4,r,2) + y**2*SorH(5,r,2))))/r**2, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: (y*z*((2 - 8*x**2)*SorH(4,r) + (-2 + 8*x**2)*SorH(5,r) + r*((-1 + 5*x**2)*SorH(4,r,1) + (1 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: (x*y*((6 - 8*x**2)*SorH(4,r) + (-6 + 8*x**2)*SorH(5,r) + r*((-3 + 5*x**2)*SorH(4,r,1) + (3 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,0,0) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*z*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x**2*z*SorH(3,r,2), \
	(1,0,1,-1) : lambda r,x,y,z,SorH: (y*z*((2 - 8*x**2)*SorH(4,r) + (-2 + 8*x**2)*SorH(5,r) + r*((-1 + 5*x**2)*SorH(4,r,1) + (1 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (2*(-1 + 4*x**2)*(-1 + x**2 + y**2)*SorH(4,r) + 2*(-1 + 4*x**2)*z**2*SorH(5,r) + r*((y**2 - 5*x**2*(-1 + x**2 + y**2))*SorH(4,r,1) + (1 - 5*x**2)*z**2*SorH(5,r,1) + r*x**2*((x**2 + y**2)*SorH(4,r,2) + z**2*SorH(5,r,2))))/r**2, \
	(1,0,1,1) : lambda r,x,y,z,SorH: (x*z*((6 - 8*x**2)*SorH(4,r) + (-6 + 8*x**2)*SorH(5,r) + r*((-3 + 5*x**2)*SorH(4,r,1) + (3 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,0,0) : lambda r,x,y,z,SorH: (3*x*(-1 + x**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x**3*SorH(3,r,2), \
	(1,1,1,-1) : lambda r,x,y,z,SorH: (x*y*((6 - 8*x**2)*SorH(4,r) + (-6 + 8*x**2)*SorH(5,r) + r*((-3 + 5*x**2)*SorH(4,r,1) + (3 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,0) : lambda r,x,y,z,SorH: (x*z*((6 - 8*x**2)*SorH(4,r) + (-6 + 8*x**2)*SorH(5,r) + r*((-3 + 5*x**2)*SorH(4,r,1) + (3 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (2*(-1 + 4*x**2)*(y**2 + z**2)*SorH(4,r) + 2*(1 - 5*x**2 + 4*x**4)*SorH(5,r) + r*(-((-1 + 5*x**2)*(y**2 + z**2)*SorH(4,r,1)) + x**2*(r*(y**2 + z**2)*SorH(4,r,2) - 5*(-1 + x**2)*SorH(5,r,1) + r*x**2*SorH(5,r,2))))/r**2, \
},
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: x*y*(-(SorH(0,r,1)/r) + SorH(0,r,2)), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: (x*(-1 + 3*y**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x*y**2*SorH(2,r,2), \
	(0,0,1,0) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(2,r) + r*(-3*SorH(2,r,1) + r*SorH(2,r,2))))/r**2, \
	(0,0,1,1) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*y*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x**2*y*SorH(2,r,2), \
	(1,-1,0,0) : lambda r,x,y,z,SorH: (x*(-1 + 3*y**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x*y**2*SorH(3,r,2), \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (x*y*((-4 + 8*x**2 + 8*z**2)*SorH(4,r) + (-4 + 8*y**2)*SorH(5,r) + r*((2 - 5*x**2 - 5*z**2)*SorH(4,r,1) + r*(x**2 + z**2)*SorH(4,r,2) + (2 - 5*y**2)*SorH(5,r,1) + r*y**2*SorH(5,r,2))))/r**2, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: (x*z*((2 - 8*y**2)*SorH(4,r) + (-2 + 8*y**2)*SorH(5,r) + r*((-1 + 5*y**2)*SorH(4,r,1) + (1 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: ((-1 + 2*y**2 + x**2*(2 - 8*y**2))*SorH(4,r) + (1 - 2*y**2 + x**2*(-2 + 8*y**2))*SorH(5,r) + r*((-y**2 + x**2*(-1 + 5*y**2))*SorH(4,r,1) + (x**2 + y**2 - 5*x**2*y**2)*SorH(5,r,1) + r*x**2*y**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,0,0,0) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(3,r) + r*(-3*SorH(3,r,1) + r*SorH(3,r,2))))/r**2, \
	(1,0,1,-1) : lambda r,x,y,z,SorH: (x*z*((2 - 8*y**2)*SorH(4,r) + (-2 + 8*y**2)*SorH(5,r) + r*((-1 + 5*y**2)*SorH(4,r,1) + (1 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (x*y*(8*(-1 + x**2 + y**2)*SorH(4,r) + 8*z**2*SorH(5,r) + r*((4 - 5*x**2 - 5*y**2)*SorH(4,r,1) + r*(x**2 + y**2)*SorH(4,r,2) + z**2*(-5*SorH(5,r,1) + r*SorH(5,r,2)))))/r**2, \
	(1,0,1,1) : lambda r,x,y,z,SorH: (y*z*((2 - 8*x**2)*SorH(4,r) + (-2 + 8*x**2)*SorH(5,r) + r*((-1 + 5*x**2)*SorH(4,r,1) + (1 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,0,0) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*y*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x**2*y*SorH(3,r,2), \
	(1,1,1,-1) : lambda r,x,y,z,SorH: ((-1 + 2*y**2 + x**2*(2 - 8*y**2))*SorH(4,r) + (1 - 2*y**2 + x**2*(-2 + 8*y**2))*SorH(5,r) + r*((-y**2 + x**2*(-1 + 5*y**2))*SorH(4,r,1) + (x**2 + y**2 - 5*x**2*y**2)*SorH(5,r,1) + r*x**2*y**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,1,1,0) : lambda r,x,y,z,SorH: (y*z*((2 - 8*x**2)*SorH(4,r) + (-2 + 8*x**2)*SorH(5,r) + r*((-1 + 5*x**2)*SorH(4,r,1) + (1 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (x*y*((-4 + 8*y**2 + 8*z**2)*SorH(4,r) + (-4 + 8*x**2)*SorH(5,r) + r*((2 - 5*y**2 - 5*z**2)*SorH(4,r,1) + r*(y**2 + z**2)*SorH(4,r,2) + (2 - 5*x**2)*SorH(5,r,1) + r*x**2*SorH(5,r,2))))/r**2, \
},
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: x*z*(-(SorH(0,r,1)/r) + SorH(0,r,2)), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(2,r) + r*(-3*SorH(2,r,1) + r*SorH(2,r,2))))/r**2, \
	(0,0,1,0) : lambda r,x,y,z,SorH: (x*(-1 + 3*z**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x*z**2*SorH(2,r,2), \
	(0,0,1,1) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*z*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x**2*z*SorH(2,r,2), \
	(1,-1,0,0) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(3,r) + r*(-3*SorH(3,r,1) + r*SorH(3,r,2))))/r**2, \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (x*z*(8*(-1 + x**2 + z**2)*SorH(4,r) + 8*y**2*SorH(5,r) + r*((4 - 5*x**2 - 5*z**2)*SorH(4,r,1) + r*(x**2 + z**2)*SorH(4,r,2) + y**2*(-5*SorH(5,r,1) + r*SorH(5,r,2)))))/r**2, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: (x*y*((2 - 8*z**2)*SorH(4,r) + (-2 + 8*z**2)*SorH(5,r) + r*((-1 + 5*z**2)*SorH(4,r,1) + (1 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: (y*z*((2 - 8*x**2)*SorH(4,r) + (-2 + 8*x**2)*SorH(5,r) + r*((-1 + 5*x**2)*SorH(4,r,1) + (1 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,0,0) : lambda r,x,y,z,SorH: (x*(-1 + 3*z**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x*z**2*SorH(3,r,2), \
	(1,0,1,-1) : lambda r,x,y,z,SorH: (x*y*((2 - 8*z**2)*SorH(4,r) + (-2 + 8*z**2)*SorH(5,r) + r*((-1 + 5*z**2)*SorH(4,r,1) + (1 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (x*z*((-4 + 8*x**2 + 8*y**2)*SorH(4,r) + (-4 + 8*z**2)*SorH(5,r) + r*((2 - 5*x**2 - 5*y**2)*SorH(4,r,1) + r*(x**2 + y**2)*SorH(4,r,2) + (2 - 5*z**2)*SorH(5,r,1) + r*z**2*SorH(5,r,2))))/r**2, \
	(1,0,1,1) : lambda r,x,y,z,SorH: ((-1 + 2*z**2 + x**2*(2 - 8*z**2))*SorH(4,r) + (1 - 2*z**2 + x**2*(-2 + 8*z**2))*SorH(5,r) + r*((-z**2 + x**2*(-1 + 5*z**2))*SorH(4,r,1) + (x**2 + z**2 - 5*x**2*z**2)*SorH(5,r,1) + r*x**2*z**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,1,0,0) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*z*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x**2*z*SorH(3,r,2), \
	(1,1,1,-1) : lambda r,x,y,z,SorH: (y*z*((2 - 8*x**2)*SorH(4,r) + (-2 + 8*x**2)*SorH(5,r) + r*((-1 + 5*x**2)*SorH(4,r,1) + (1 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,0) : lambda r,x,y,z,SorH: ((-1 + 2*z**2 + x**2*(2 - 8*z**2))*SorH(4,r) + (1 - 2*z**2 + x**2*(-2 + 8*z**2))*SorH(5,r) + r*((-z**2 + x**2*(-1 + 5*z**2))*SorH(4,r,1) + (x**2 + z**2 - 5*x**2*z**2)*SorH(5,r,1) + r*x**2*z**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (x*z*((-4 + 8*y**2 + 8*z**2)*SorH(4,r) + (-4 + 8*x**2)*SorH(5,r) + r*((2 - 5*y**2 - 5*z**2)*SorH(4,r,1) + r*(y**2 + z**2)*SorH(4,r,2) + (2 - 5*x**2)*SorH(5,r,1) + r*x**2*SorH(5,r,2))))/r**2, \
},
],\
[\
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: x*y*(-(SorH(0,r,1)/r) + SorH(0,r,2)), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: (x*(-1 + 3*y**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x*y**2*SorH(2,r,2), \
	(0,0,1,0) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(2,r) + r*(-3*SorH(2,r,1) + r*SorH(2,r,2))))/r**2, \
	(0,0,1,1) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*y*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x**2*y*SorH(2,r,2), \
	(1,-1,0,0) : lambda r,x,y,z,SorH: (x*(-1 + 3*y**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x*y**2*SorH(3,r,2), \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (x*y*((-4 + 8*x**2 + 8*z**2)*SorH(4,r) + (-4 + 8*y**2)*SorH(5,r) + r*((2 - 5*x**2 - 5*z**2)*SorH(4,r,1) + r*(x**2 + z**2)*SorH(4,r,2) + (2 - 5*y**2)*SorH(5,r,1) + r*y**2*SorH(5,r,2))))/r**2, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: (x*z*((2 - 8*y**2)*SorH(4,r) + (-2 + 8*y**2)*SorH(5,r) + r*((-1 + 5*y**2)*SorH(4,r,1) + (1 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: ((-1 + 2*y**2 + x**2*(2 - 8*y**2))*SorH(4,r) + (1 - 2*y**2 + x**2*(-2 + 8*y**2))*SorH(5,r) + r*((-y**2 + x**2*(-1 + 5*y**2))*SorH(4,r,1) + (x**2 + y**2 - 5*x**2*y**2)*SorH(5,r,1) + r*x**2*y**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,0,0,0) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(3,r) + r*(-3*SorH(3,r,1) + r*SorH(3,r,2))))/r**2, \
	(1,0,1,-1) : lambda r,x,y,z,SorH: (x*z*((2 - 8*y**2)*SorH(4,r) + (-2 + 8*y**2)*SorH(5,r) + r*((-1 + 5*y**2)*SorH(4,r,1) + (1 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (x*y*(8*(-1 + x**2 + y**2)*SorH(4,r) + 8*z**2*SorH(5,r) + r*((4 - 5*x**2 - 5*y**2)*SorH(4,r,1) + r*(x**2 + y**2)*SorH(4,r,2) + z**2*(-5*SorH(5,r,1) + r*SorH(5,r,2)))))/r**2, \
	(1,0,1,1) : lambda r,x,y,z,SorH: (y*z*((2 - 8*x**2)*SorH(4,r) + (-2 + 8*x**2)*SorH(5,r) + r*((-1 + 5*x**2)*SorH(4,r,1) + (1 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,0,0) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*y*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x**2*y*SorH(3,r,2), \
	(1,1,1,-1) : lambda r,x,y,z,SorH: ((-1 + 2*y**2 + x**2*(2 - 8*y**2))*SorH(4,r) + (1 - 2*y**2 + x**2*(-2 + 8*y**2))*SorH(5,r) + r*((-y**2 + x**2*(-1 + 5*y**2))*SorH(4,r,1) + (x**2 + y**2 - 5*x**2*y**2)*SorH(5,r,1) + r*x**2*y**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,1,1,0) : lambda r,x,y,z,SorH: (y*z*((2 - 8*x**2)*SorH(4,r) + (-2 + 8*x**2)*SorH(5,r) + r*((-1 + 5*x**2)*SorH(4,r,1) + (1 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (x*y*((-4 + 8*y**2 + 8*z**2)*SorH(4,r) + (-4 + 8*x**2)*SorH(5,r) + r*((2 - 5*y**2 - 5*z**2)*SorH(4,r,1) + r*(y**2 + z**2)*SorH(4,r,2) + (2 - 5*x**2)*SorH(5,r,1) + r*x**2*SorH(5,r,2))))/r**2, \
},
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: -(((-1 + y**2)*SorH(0,r,1))/r) + y**2*SorH(0,r,2), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: (3*y*(-1 + y**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + y**3*SorH(2,r,2), \
	(0,0,1,0) : lambda r,x,y,z,SorH: ((-1 + 3*y**2)*z*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + y**2*z*SorH(2,r,2), \
	(0,0,1,1) : lambda r,x,y,z,SorH: (x*(-1 + 3*y**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x*y**2*SorH(2,r,2), \
	(1,-1,0,0) : lambda r,x,y,z,SorH: (3*y*(-1 + y**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + y**3*SorH(3,r,2), \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (2*(-1 + 4*y**2)*(x**2 + z**2)*SorH(4,r) + 2*(1 - 5*y**2 + 4*y**4)*SorH(5,r) + r*(-((-1 + 5*y**2)*(x**2 + z**2)*SorH(4,r,1)) + y**2*(r*(x**2 + z**2)*SorH(4,r,2) - 5*(-1 + y**2)*SorH(5,r,1) + r*y**2*SorH(5,r,2))))/r**2, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: (y*z*((6 - 8*y**2)*SorH(4,r) + (-6 + 8*y**2)*SorH(5,r) + r*((-3 + 5*y**2)*SorH(4,r,1) + (3 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: (x*y*((6 - 8*y**2)*SorH(4,r) + (-6 + 8*y**2)*SorH(5,r) + r*((-3 + 5*y**2)*SorH(4,r,1) + (3 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,0,0) : lambda r,x,y,z,SorH: ((-1 + 3*y**2)*z*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + y**2*z*SorH(3,r,2), \
	(1,0,1,-1) : lambda r,x,y,z,SorH: (y*z*((6 - 8*y**2)*SorH(4,r) + (-6 + 8*y**2)*SorH(5,r) + r*((-3 + 5*y**2)*SorH(4,r,1) + (3 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (2*(-1 + x**2 + y**2)*(-1 + 4*y**2)*SorH(4,r) + 2*(-1 + 4*y**2)*z**2*SorH(5,r) + r*((x**2 - 5*(-1 + x**2)*y**2 - 5*y**4)*SorH(4,r,1) + (1 - 5*y**2)*z**2*SorH(5,r,1) + r*y**2*((x**2 + y**2)*SorH(4,r,2) + z**2*SorH(5,r,2))))/r**2, \
	(1,0,1,1) : lambda r,x,y,z,SorH: (x*z*((2 - 8*y**2)*SorH(4,r) + (-2 + 8*y**2)*SorH(5,r) + r*((-1 + 5*y**2)*SorH(4,r,1) + (1 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,0,0) : lambda r,x,y,z,SorH: (x*(-1 + 3*y**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x*y**2*SorH(3,r,2), \
	(1,1,1,-1) : lambda r,x,y,z,SorH: (x*y*((6 - 8*y**2)*SorH(4,r) + (-6 + 8*y**2)*SorH(5,r) + r*((-3 + 5*y**2)*SorH(4,r,1) + (3 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,0) : lambda r,x,y,z,SorH: (x*z*((2 - 8*y**2)*SorH(4,r) + (-2 + 8*y**2)*SorH(5,r) + r*((-1 + 5*y**2)*SorH(4,r,1) + (1 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (2*(-1 + 4*y**2)*(-1 + y**2 + z**2)*SorH(4,r) + 2*x**2*(-1 + 4*y**2)*SorH(5,r) + r*((z**2 - 5*y**2*(-1 + y**2 + z**2))*SorH(4,r,1) + x**2*(1 - 5*y**2)*SorH(5,r,1) + r*y**2*((y**2 + z**2)*SorH(4,r,2) + x**2*SorH(5,r,2))))/r**2, \
},
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: y*z*(-(SorH(0,r,1)/r) + SorH(0,r,2)), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: ((-1 + 3*y**2)*z*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + y**2*z*SorH(2,r,2), \
	(0,0,1,0) : lambda r,x,y,z,SorH: (y*(-1 + 3*z**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + y*z**2*SorH(2,r,2), \
	(0,0,1,1) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(2,r) + r*(-3*SorH(2,r,1) + r*SorH(2,r,2))))/r**2, \
	(1,-1,0,0) : lambda r,x,y,z,SorH: ((-1 + 3*y**2)*z*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + y**2*z*SorH(3,r,2), \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (y*z*((-4 + 8*x**2 + 8*z**2)*SorH(4,r) + (-4 + 8*y**2)*SorH(5,r) + r*((2 - 5*x**2 - 5*z**2)*SorH(4,r,1) + r*(x**2 + z**2)*SorH(4,r,2) + (2 - 5*y**2)*SorH(5,r,1) + r*y**2*SorH(5,r,2))))/r**2, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: ((-1 + 2*z**2 + y**2*(2 - 8*z**2))*SorH(4,r) + (1 - 2*z**2 + y**2*(-2 + 8*z**2))*SorH(5,r) + r*((-z**2 + y**2*(-1 + 5*z**2))*SorH(4,r,1) + (y**2 + z**2 - 5*y**2*z**2)*SorH(5,r,1) + r*y**2*z**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: (x*z*((2 - 8*y**2)*SorH(4,r) + (-2 + 8*y**2)*SorH(5,r) + r*((-1 + 5*y**2)*SorH(4,r,1) + (1 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,0,0) : lambda r,x,y,z,SorH: (y*(-1 + 3*z**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + y*z**2*SorH(3,r,2), \
	(1,0,1,-1) : lambda r,x,y,z,SorH: ((-1 + 2*z**2 + y**2*(2 - 8*z**2))*SorH(4,r) + (1 - 2*z**2 + y**2*(-2 + 8*z**2))*SorH(5,r) + r*((-z**2 + y**2*(-1 + 5*z**2))*SorH(4,r,1) + (y**2 + z**2 - 5*y**2*z**2)*SorH(5,r,1) + r*y**2*z**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (y*z*((-4 + 8*x**2 + 8*y**2)*SorH(4,r) + (-4 + 8*z**2)*SorH(5,r) + r*((2 - 5*x**2 - 5*y**2)*SorH(4,r,1) + r*(x**2 + y**2)*SorH(4,r,2) + (2 - 5*z**2)*SorH(5,r,1) + r*z**2*SorH(5,r,2))))/r**2, \
	(1,0,1,1) : lambda r,x,y,z,SorH: (x*y*((2 - 8*z**2)*SorH(4,r) + (-2 + 8*z**2)*SorH(5,r) + r*((-1 + 5*z**2)*SorH(4,r,1) + (1 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,0,0) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(3,r) + r*(-3*SorH(3,r,1) + r*SorH(3,r,2))))/r**2, \
	(1,1,1,-1) : lambda r,x,y,z,SorH: (x*z*((2 - 8*y**2)*SorH(4,r) + (-2 + 8*y**2)*SorH(5,r) + r*((-1 + 5*y**2)*SorH(4,r,1) + (1 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,0) : lambda r,x,y,z,SorH: (x*y*((2 - 8*z**2)*SorH(4,r) + (-2 + 8*z**2)*SorH(5,r) + r*((-1 + 5*z**2)*SorH(4,r,1) + (1 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (y*z*(8*(-1 + y**2 + z**2)*SorH(4,r) + 8*x**2*SorH(5,r) + r*((4 - 5*y**2 - 5*z**2)*SorH(4,r,1) + r*(y**2 + z**2)*SorH(4,r,2) + x**2*(-5*SorH(5,r,1) + r*SorH(5,r,2)))))/r**2, \
},
],\
[\
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: x*z*(-(SorH(0,r,1)/r) + SorH(0,r,2)), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(2,r) + r*(-3*SorH(2,r,1) + r*SorH(2,r,2))))/r**2, \
	(0,0,1,0) : lambda r,x,y,z,SorH: (x*(-1 + 3*z**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x*z**2*SorH(2,r,2), \
	(0,0,1,1) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*z*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x**2*z*SorH(2,r,2), \
	(1,-1,0,0) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(3,r) + r*(-3*SorH(3,r,1) + r*SorH(3,r,2))))/r**2, \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (x*z*(8*(-1 + x**2 + z**2)*SorH(4,r) + 8*y**2*SorH(5,r) + r*((4 - 5*x**2 - 5*z**2)*SorH(4,r,1) + r*(x**2 + z**2)*SorH(4,r,2) + y**2*(-5*SorH(5,r,1) + r*SorH(5,r,2)))))/r**2, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: (x*y*((2 - 8*z**2)*SorH(4,r) + (-2 + 8*z**2)*SorH(5,r) + r*((-1 + 5*z**2)*SorH(4,r,1) + (1 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: (y*z*((2 - 8*x**2)*SorH(4,r) + (-2 + 8*x**2)*SorH(5,r) + r*((-1 + 5*x**2)*SorH(4,r,1) + (1 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,0,0) : lambda r,x,y,z,SorH: (x*(-1 + 3*z**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x*z**2*SorH(3,r,2), \
	(1,0,1,-1) : lambda r,x,y,z,SorH: (x*y*((2 - 8*z**2)*SorH(4,r) + (-2 + 8*z**2)*SorH(5,r) + r*((-1 + 5*z**2)*SorH(4,r,1) + (1 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (x*z*((-4 + 8*x**2 + 8*y**2)*SorH(4,r) + (-4 + 8*z**2)*SorH(5,r) + r*((2 - 5*x**2 - 5*y**2)*SorH(4,r,1) + r*(x**2 + y**2)*SorH(4,r,2) + (2 - 5*z**2)*SorH(5,r,1) + r*z**2*SorH(5,r,2))))/r**2, \
	(1,0,1,1) : lambda r,x,y,z,SorH: ((-1 + 2*z**2 + x**2*(2 - 8*z**2))*SorH(4,r) + (1 - 2*z**2 + x**2*(-2 + 8*z**2))*SorH(5,r) + r*((-z**2 + x**2*(-1 + 5*z**2))*SorH(4,r,1) + (x**2 + z**2 - 5*x**2*z**2)*SorH(5,r,1) + r*x**2*z**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,1,0,0) : lambda r,x,y,z,SorH: ((-1 + 3*x**2)*z*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x**2*z*SorH(3,r,2), \
	(1,1,1,-1) : lambda r,x,y,z,SorH: (y*z*((2 - 8*x**2)*SorH(4,r) + (-2 + 8*x**2)*SorH(5,r) + r*((-1 + 5*x**2)*SorH(4,r,1) + (1 - 5*x**2)*SorH(5,r,1) + r*x**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,0) : lambda r,x,y,z,SorH: ((-1 + 2*z**2 + x**2*(2 - 8*z**2))*SorH(4,r) + (1 - 2*z**2 + x**2*(-2 + 8*z**2))*SorH(5,r) + r*((-z**2 + x**2*(-1 + 5*z**2))*SorH(4,r,1) + (x**2 + z**2 - 5*x**2*z**2)*SorH(5,r,1) + r*x**2*z**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (x*z*((-4 + 8*y**2 + 8*z**2)*SorH(4,r) + (-4 + 8*x**2)*SorH(5,r) + r*((2 - 5*y**2 - 5*z**2)*SorH(4,r,1) + r*(y**2 + z**2)*SorH(4,r,2) + (2 - 5*x**2)*SorH(5,r,1) + r*x**2*SorH(5,r,2))))/r**2, \
},
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: y*z*(-(SorH(0,r,1)/r) + SorH(0,r,2)), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: ((-1 + 3*y**2)*z*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + y**2*z*SorH(2,r,2), \
	(0,0,1,0) : lambda r,x,y,z,SorH: (y*(-1 + 3*z**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + y*z**2*SorH(2,r,2), \
	(0,0,1,1) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(2,r) + r*(-3*SorH(2,r,1) + r*SorH(2,r,2))))/r**2, \
	(1,-1,0,0) : lambda r,x,y,z,SorH: ((-1 + 3*y**2)*z*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + y**2*z*SorH(3,r,2), \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (y*z*((-4 + 8*x**2 + 8*z**2)*SorH(4,r) + (-4 + 8*y**2)*SorH(5,r) + r*((2 - 5*x**2 - 5*z**2)*SorH(4,r,1) + r*(x**2 + z**2)*SorH(4,r,2) + (2 - 5*y**2)*SorH(5,r,1) + r*y**2*SorH(5,r,2))))/r**2, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: ((-1 + 2*z**2 + y**2*(2 - 8*z**2))*SorH(4,r) + (1 - 2*z**2 + y**2*(-2 + 8*z**2))*SorH(5,r) + r*((-z**2 + y**2*(-1 + 5*z**2))*SorH(4,r,1) + (y**2 + z**2 - 5*y**2*z**2)*SorH(5,r,1) + r*y**2*z**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: (x*z*((2 - 8*y**2)*SorH(4,r) + (-2 + 8*y**2)*SorH(5,r) + r*((-1 + 5*y**2)*SorH(4,r,1) + (1 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,0,0) : lambda r,x,y,z,SorH: (y*(-1 + 3*z**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + y*z**2*SorH(3,r,2), \
	(1,0,1,-1) : lambda r,x,y,z,SorH: ((-1 + 2*z**2 + y**2*(2 - 8*z**2))*SorH(4,r) + (1 - 2*z**2 + y**2*(-2 + 8*z**2))*SorH(5,r) + r*((-z**2 + y**2*(-1 + 5*z**2))*SorH(4,r,1) + (y**2 + z**2 - 5*y**2*z**2)*SorH(5,r,1) + r*y**2*z**2*(-SorH(4,r,2) + SorH(5,r,2))))/r**2, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (y*z*((-4 + 8*x**2 + 8*y**2)*SorH(4,r) + (-4 + 8*z**2)*SorH(5,r) + r*((2 - 5*x**2 - 5*y**2)*SorH(4,r,1) + r*(x**2 + y**2)*SorH(4,r,2) + (2 - 5*z**2)*SorH(5,r,1) + r*z**2*SorH(5,r,2))))/r**2, \
	(1,0,1,1) : lambda r,x,y,z,SorH: (x*y*((2 - 8*z**2)*SorH(4,r) + (-2 + 8*z**2)*SorH(5,r) + r*((-1 + 5*z**2)*SorH(4,r,1) + (1 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,0,0) : lambda r,x,y,z,SorH: (x*y*z*(3*SorH(3,r) + r*(-3*SorH(3,r,1) + r*SorH(3,r,2))))/r**2, \
	(1,1,1,-1) : lambda r,x,y,z,SorH: (x*z*((2 - 8*y**2)*SorH(4,r) + (-2 + 8*y**2)*SorH(5,r) + r*((-1 + 5*y**2)*SorH(4,r,1) + (1 - 5*y**2)*SorH(5,r,1) + r*y**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,0) : lambda r,x,y,z,SorH: (x*y*((2 - 8*z**2)*SorH(4,r) + (-2 + 8*z**2)*SorH(5,r) + r*((-1 + 5*z**2)*SorH(4,r,1) + (1 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (y*z*(8*(-1 + y**2 + z**2)*SorH(4,r) + 8*x**2*SorH(5,r) + r*((4 - 5*y**2 - 5*z**2)*SorH(4,r,1) + r*(y**2 + z**2)*SorH(4,r,2) + x**2*(-5*SorH(5,r,1) + r*SorH(5,r,2)))))/r**2, \
},
{\
	(0,0,0,0) : lambda r,x,y,z,SorH: -(((-1 + z**2)*SorH(0,r,1))/r) + z**2*SorH(0,r,2), \
	(0,0,1,-1) : lambda r,x,y,z,SorH: (y*(-1 + 3*z**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + y*z**2*SorH(2,r,2), \
	(0,0,1,0) : lambda r,x,y,z,SorH: (3*z*(-1 + z**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + z**3*SorH(2,r,2), \
	(0,0,1,1) : lambda r,x,y,z,SorH: (x*(-1 + 3*z**2)*(SorH(2,r) - r*SorH(2,r,1)))/r**2 + x*z**2*SorH(2,r,2), \
	(1,-1,0,0) : lambda r,x,y,z,SorH: (y*(-1 + 3*z**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + y*z**2*SorH(3,r,2), \
	(1,-1,1,-1) : lambda r,x,y,z,SorH: (2*(-1 + x**2 + z**2)*(-1 + 4*z**2)*SorH(4,r) + 2*y**2*(-1 + 4*z**2)*SorH(5,r) + r*((x**2 - 5*(-1 + x**2)*z**2 - 5*z**4)*SorH(4,r,1) + y**2*(1 - 5*z**2)*SorH(5,r,1) + r*z**2*((x**2 + z**2)*SorH(4,r,2) + y**2*SorH(5,r,2))))/r**2, \
	(1,-1,1,0) : lambda r,x,y,z,SorH: (y*z*((6 - 8*z**2)*SorH(4,r) + (-6 + 8*z**2)*SorH(5,r) + r*((-3 + 5*z**2)*SorH(4,r,1) + (3 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,-1,1,1) : lambda r,x,y,z,SorH: (x*y*((2 - 8*z**2)*SorH(4,r) + (-2 + 8*z**2)*SorH(5,r) + r*((-1 + 5*z**2)*SorH(4,r,1) + (1 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,0,0) : lambda r,x,y,z,SorH: (3*z*(-1 + z**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + z**3*SorH(3,r,2), \
	(1,0,1,-1) : lambda r,x,y,z,SorH: (y*z*((6 - 8*z**2)*SorH(4,r) + (-6 + 8*z**2)*SorH(5,r) + r*((-3 + 5*z**2)*SorH(4,r,1) + (3 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,0,1,0) : lambda r,x,y,z,SorH: (2*(x**2 + y**2)*(-1 + 4*z**2)*SorH(4,r) + 2*(1 - 5*z**2 + 4*z**4)*SorH(5,r) + r*(-((x**2 + y**2)*(-1 + 5*z**2)*SorH(4,r,1)) + z**2*(r*(x**2 + y**2)*SorH(4,r,2) - 5*(-1 + z**2)*SorH(5,r,1) + r*z**2*SorH(5,r,2))))/r**2, \
	(1,0,1,1) : lambda r,x,y,z,SorH: (x*z*((6 - 8*z**2)*SorH(4,r) + (-6 + 8*z**2)*SorH(5,r) + r*((-3 + 5*z**2)*SorH(4,r,1) + (3 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,0,0) : lambda r,x,y,z,SorH: (x*(-1 + 3*z**2)*(SorH(3,r) - r*SorH(3,r,1)))/r**2 + x*z**2*SorH(3,r,2), \
	(1,1,1,-1) : lambda r,x,y,z,SorH: (x*y*((2 - 8*z**2)*SorH(4,r) + (-2 + 8*z**2)*SorH(5,r) + r*((-1 + 5*z**2)*SorH(4,r,1) + (1 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,0) : lambda r,x,y,z,SorH: (x*z*((6 - 8*z**2)*SorH(4,r) + (-6 + 8*z**2)*SorH(5,r) + r*((-3 + 5*z**2)*SorH(4,r,1) + (3 - 5*z**2)*SorH(5,r,1) + r*z**2*(-SorH(4,r,2) + SorH(5,r,2)))))/r**2, \
	(1,1,1,1) : lambda r,x,y,z,SorH: (2*(-1 + y**2 + z**2)*(-1 + 4*z**2)*SorH(4,r) + 2*x**2*(-1 + 4*z**2)*SorH(5,r) + r*((y**2 - 5*(-1 + y**2)*z**2 - 5*z**4)*SorH(4,r,1) + x**2*(1 - 5*z**2)*SorH(5,r,1) + r*z**2*((y**2 + z**2)*SorH(4,r,2) + x**2*SorH(5,r,2))))/r**2, \
},
],\
]
