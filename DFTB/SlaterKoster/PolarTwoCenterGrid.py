"""
Routines for integrating functions that have peaks at two different centers on 
a double polar grid.
"""

from numpy import linspace, pi, sin, cos, tan, arccos, exp, log, sqrt, array, inf, ediff1d, tile, vstack
from itertools import product

def ptcgrid(hs,rs,angles):
    """
    Generate a polar two center grid for each separation of the centers.

    Parameters:
    ===========
    hs: 1D numpy array, locations of centers at y=-h and y=+h
    rs: 1D numpy array, radial grid covering [0, rmax]
    angles: 1D numpy array, angular grid covering interval [0,pi]

    Returns:
    ========
    Xs: list of 1D numpy arrays, Xs[i] contains the x-positions of the grid 
        points for a grid where the centers are separated by a distance 2*hs[i]
    Ys: list of 1D numpy arrays, Ys[i] contains the y-positions for a separation 2*hs[i]
    Areas: list of 1D numpy arrays, Areas[i] contains the areas of the quads for a separation 2*hs[i]

    Usage:
    ======
    You can find an integral as a function of the separation of the two centers by
      I = [sum(f(Xs[i],Ys[i],h)*Areas[i]) for i,h in enumerate(hs)]
    """
    assert rs[0] > 0.0 # cannot cover the very center
    Xs,Ys,Areas = [], [], []
    for h in hs:
        print "constructing grid for h=%s" % h
        p = PolarTwoCenterGrid(h)
        p.setGrid(rs,angles)
        p.tessellate(halfplanes="+")
        xs,ys,areas = p.getQuads()
        Xs.append(xs)
        Ys.append(ys)
        Areas.append(areas)
    return Xs, Ys, Areas
        
class PolarTwoCenterGrid:
    def __init__(self, h):
        """h is the separation between the two centers along the y-axis"""
        self.h = h
        self.centers_x = []
        self.centers_y = []
        self.areas = []
    def setGrid(self, r, angle):
        """
        Set the grid resolution and grid spacing.

        Parameters:
        ===========
        r: 1D numpy array, radial grid points which do not have to be equidistant
        angle: 1D numpy array, angular grid points covering the interval [0,pi]
        """
        assert len(angle) >= 6 # otherwise the tessalation is wrong
        assert angle[0] == 0.0 and angle[-1] == pi
        self.r = r
        self.angle = angle
        self.dr = ediff1d(self.r, to_end=self.r[-1]-self.r[-2])
        self.da = ediff1d(self.angle, to_end=self.angle[-1]-self.angle[-2])
    def tessellate(self, halfplanes="+-"):
        """
        compute centers and areas of the grid and store them internally. Retrieve
        the grid points with the method getQuads().

        Parameters:
        ===========
        halfplane (optional, default "+-"): 
           "+" compute only those grid points that lie in the x>=0 half plane
               (for integrands that are spherically symmetric about y-axis)
           "+-" compute grid points in both half planes
        """
        self.halfplanes = halfplanes
        def polar2cart(r,a):
            """polar coordinates relative to a center that is shifted
            up along the z-axis by h."""
            x = r*sin(a)
            y = r*cos(a) + self.h
            return x,y
        h = self.h
        for (ri,ai),(dr,da) in zip(product(self.r,self.angle), product(self.dr, self.da)):
            xi,yi = polar2cart(ri,ai)
            if yi < 0.0:
                # only upper quadrant
                continue
            if ai == pi:
                # do not duplicate quads that lie on the y-axis
                continue
            corners = []
            def RR0(a):
                """find the radial component of the line y=0
                in polar coordinates for angle a"""
                if a <= pi/2.0:
                    # for a < pi/2 the radial vector never intersects y=0
                    return inf
                return abs(h/cos(a))
            def AA0(r):
                """find angle for which the radial vector of length r
                intersects the line y=0"""
                if ai <= pi/2.0 or r<=h:
                    # hack
                    return pi
                elif pi/2.0 < ai <= pi:
                    a0 = arccos(-h/r)
                return a0
            # corners in shifted polar coordinates
            quad_type = ""
            corners.append((ri,ai)) # A
            quad_type += "A"
            corners.append((min(ri+dr, RR0(ai)), ai)) # B
            quad_type += "B"
            if ai <= pi/2.0 or ai < AA0(ri+dr):
                corners.append((ri+dr, min(ai+da, max(ai,AA0(ri+dr))))) # C
                quad_type += "C"
            if ri < RR0(ai+da) < ri+dr:
                corners.append((max(ri, min(ri+dr, RR0(ai+da))), ai+da)) # D
                quad_type += "D"
            corners.append((ri, min(AA0(ri), ai+da))) # E
            quad_type += "E"
            #print quad_type
            polygon_x = []
            polygon_y = []
            for (rc,ac) in corners:
                x,y = polar2cart(rc,ac)
                polygon_x.append(x)
                polygon_y.append(y)
            # center position is average of corners
            center_x = sum(polygon_x)/len(polygon_x)
            center_y = sum(polygon_y)/len(polygon_y)
            # areas
            if ri+dr < h:
                area = 0.5*da*(2.0*ri*dr + dr**2)
            else:
                acos1 = arccos(-min(h/ri,1.0))
                acos2 = arccos(-h/(ri+dr))
                tan1 = tan(acos1)
                tan2 = tan(acos2)

                if quad_type == "ABE":
                    area = h**2/2.0 * (tan1 - tan(ai)) - 0.5*ri**2 * (acos1 - ai)
                elif quad_type == "ABCE":
                    if RR0(ai+da) > ri+dr:
                        area = 0.5*da*(2.0*ri*dr + dr**2)
                    else:
                        area1 = 0.5*(2*ri*dr + dr**2)*(acos2 - ai)
                        area2 = h**2/2.0 * (tan1 - tan2) - 0.5*ri**2 * (acos1 - acos2)
                        area = area1 + area2
                elif quad_type == "ABCDE":
                    area1 = 0.5*(2*ri*dr + dr**2)*(acos2 - ai)
                    area2 = h**2/2.0 * (tan(ai+da) - tan2) - 0.5*ri**2 * (ai+da - acos2)
                    area = area1 + area2
                elif quad_type == "ABDE":
                    assert (RR0(ai) < ri+dr and RR0(ai+da) < ri+dr)
                    area = h**2/2.0 * (tan(ai+da) - tan(ai)) - 0.5*ri**2 * da
                else:
                    raise Exception("???")
            """
            # uncomment this block to plot the grid
            color_code = {"ABE": "red", "ABCE": "black", "ABDE": "blue", "ABCDE": "green"}
            # close polygon
            polygon_x.append(polygon_x[0])
            polygon_y.append(polygon_y[0])
 
            plot([center_x],  [center_y], "o", color="red")
            plot([center_x],  [-center_y], "o", color="red")
            plot([-center_x], [-center_y], "o", color="red")
            plot([-center_x], [center_y], "o", color="red")

            plot(polygon_x, polygon_y, color=color_code[quad_type])
#            if area > 0.2:
#                text(center_x, center_y, "%2.3f" % area)
            """
            # The centers of the other three quadrants are obtained by reflection
            # on the x- and y-axis of the centers in the upper right quadrant
            # upper right quadrant, lower right quadrant, lower left quadrant, upper left quadrant
            self.centers_x += [center_x,  center_x]
            self.centers_y += [center_y, -center_y]
            self.areas     += [area, area]
            if self.halfplanes == "+-":
                self.centers_x += [-center_x, -center_x]
                self.centers_y += [-center_y,  center_y]
                self.areas     += [area, area]

        self.centers_x = array(self.centers_x)
        self.centers_y = array(self.centers_y)
        self.areas = array(self.areas)
        # self test
        I1 = sum(self.areas)
        I2 = _total_grid_area(self.h, self.r[-1], self.r[-1]-self.r[-2])
        if halfplanes != "+-":
            I2 *= 0.5
        print "numerical area  = %s" % I1
        print "total grid area = %s" % I2
        #assert abs(I1 - I2) < 1.0e-10
    def quads_it(self):
        """gives iterator to centers and areas"""
        for cx,cy,da in zip(self.centers_x, self.centers_y, self.areas):
            yield cx,cy,da
    def getQuads(self):
        """
        gives centers of quads and their exact areas.
        You can integrate a function defined in the x-y-plane like this:
         I = sum(f(x,y)*areas)

        Returns:
        ========
        x: 1D numpy array, x-positions of quad centers
        y: 1D numpy array, y-positions of quad centers
        areas: 1D numpy array, exact areas of quads
        """
        return self.centers_x, self.centers_y, self.areas

def _total_grid_area(h, rmax, dr): # for testing purposes
    """
    calculate the area covered by the overlapping disks with radii (rmax+dr)
    centered at z=+h and z=-h.
    We should obtain exactly this number if we  integrate the constant function
    f(rho,z) = 1 on the two center polar grid.
    """
    R = rmax+dr
    A = 2.0*pi*R**2
    if h < R:
        a = arccos(h/R)
        A -= 2*a*R**2 - 2*h*R*sqrt(1-cos(a)**2)
    return A

__all__ = ["ptcgrid", "PolarTwoCenterGrid"]

## TESTS ##
def test_ptcgrid():
    def f(x,y,h):
        return exp(-(pow(x,2)+pow(y-h,2)))*exp(-(pow(x,2)+pow(y+h,2)))
    rmax = 10.0
    Nr = 100
    Na = 30
    angles = linspace(0.0, pi, Na)
    hs = linspace(0.0, 5.0, 200)
    
    rs = exp(linspace(log(0.00000001), log(rmax), Nr))
    Xs,Ys,Areas = ptcgrid(hs,rs,angles)
    plot(Xs[-1],Ys[-1], "o")
    show()
    I = [sum(f(Xs[i],Ys[i],h)*Areas[i]) for i,h in enumerate(hs)]
    print I
    cla()
    plot(hs,I)
    show()

def test_PolarTwoCenterGrid():
    rmax = 3.0          
    Nr = 40
    Na = 23
    h = 1.0
    p = PolarTwoCenterGrid(h)
    p.setGrid(exp(linspace(log(0.000001), log(rmax), Nr)), linspace(0.0, pi, Na))
    p.tessellate()

    def f(x,y):
        return 1.0

    xs,ys,areas = p.getQuads()
    I = sum(f(xs,ys)*areas)

    show()

if __name__ == "__main__":
    from matplotlib.pyplot import plot, cla, show, text
    test_PolarTwoCenterGrid()
#    test_ptcgrid()
