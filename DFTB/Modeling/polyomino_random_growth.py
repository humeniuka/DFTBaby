#!/usr/bin/env python
"""
perform a Monte Carlo simulation of the grow process for 
triply fused porphyrin flakes
"""

import numpy as np
import numpy.linalg as la

import numpy.random 
import matplotlib.pyplot as plt

from DFTB.Modeling.porphyrin_flakes import Flake, Gouterman_matelems

def random_selection(p):
    """
    draw integer i with probability p[i]

    Parameters
    ----------
    p           :  numpy array with probabilities for events
    
    Returns
    -------
    i           :  random index in the range 0 <= i <= len(p) with
                   probability p[i]
    """
    n = len(p)
    # normalize probabilities, sum_i p[i] = 1
    p /= np.sum(p)
    # cumulative probability distribution
    s = np.cumsum(p)
    # draw random number from interval [0,1)
    r = np.random.rand(1)[0]
    for i in range(0, n):
        if i == 0:
            a = 0.0
        else:
            a = s[i-1]
        b = s[i]
        # r lies in the interval [s[i-1], s[i]] with probability
        # p[i] sind s[i]-s[i-1] = p[i]
        assert abs((b-a) - p[i]) < 1.0e-10
        if a <= r < b:
            break
        
    return i

def grid_flood_fill_recursive(grid, i0,j0, target_color=0, replacement_color=2):
    """
    Flood fill a grid starting from the node (i,j) by replacing the colors
    of all connected squares with a different color.

    see https://en.wikipedia.org/wiki/Flood_fill (recursive implementation)

    Parameters
    ----------
    grid      :  2d numpy array, integers
    i0,j0     :  horizontal and vertical grid positions 
    
    Optional
    --------
    target_color      : integer, select grid points having grid[i,j] == target_color
    replacement_color : integer, color for replacing selected grid points 

    Returns
    -------
    grid      :  flooded grid
    """
    grid = np.copy(grid)
    nx,ny = grid.shape
    
    def flood_fill(i,j):
        if grid[i,j] != target_color:
            return
        grid[i,j] = replacement_color
        if      i+1 < nx:
            flood_fill(i+1,j)
        if 0 <= i-1:
            flood_fill(i-1,j)
        if      j+1 < ny:
            flood_fill(i,  j+1)
        if 0 <= j-1:
            flood_fill(i,  j-1)

    flood_fill(i0,j0)

    return grid

def grid_flood_fill_queue(grid, i0,j0, target_color=0, replacement_color=2):
    """
    Flood fill a grid starting from the node (i,j) by replacing the colors
    of all connected squares with a different color.

    see https://en.wikipedia.org/wiki/Flood_fill (implementation using queue)

    Parameters
    ----------
    grid      :  2d numpy array, integers
    i0,j0     :  horizontal and vertical grid positions 
    
    Optional
    --------
    target_color      : integer, select grid points having grid[i,j] == target_color
    replacement_color : integer, color for replacing selected grid points 

    Returns
    -------
    grid      :  flooded grid
    """
    grid = np.copy(grid)
    nx,ny = grid.shape

    if grid[i0,j0] != target_color:
        return

    # 
    queue = [(i0,j0)]
    while len(queue) > 0:
        i,j = queue.pop(0)
        if grid[i,j] == target_color:
            grid[i,j] = replacement_color
            # neighbouring grid points
            if      i+1 < nx:
                queue.append((i+1,j))
            if 0 <= i-1:
                queue.append((i-1,j))
            if      j+1 < ny:
                queue.append((i,  j+1))
            if 0 <= j-1:
                queue.append((i,  j-1))

    return grid

#
# Which implementation should be used for flood filling?
#
#grid_flood_fill = grid_flood_fill_recursive
grid_flood_fill = grid_flood_fill_queue


def identify_holes(grid):
    """
    identify the holes in a porphyrin flakes. These are empty sites that lie in the
    interior of an island and cannot be filled in subsequent steps, since porphyrin
    monomers would have to dive into solution to access the hole.
    
    Grid points in the hole are set to -1.
    """
    # start flood filling from upper left corner. All grid points having color
    # 0 and can be reached from the upper corner are set to 2.
    grid = grid_flood_fill(grid, 0,0, target_color=0, replacement_color=2)

    # Exterior points have color 2, points that still have color 0 must
    # be holes.
    # set color of holes to -1
    grid[grid == 0] = -1
    # revert color of exterior points back to 0
    grid[grid == 2] = 0

    return grid

def sizes_of_holes(grid):
    """
    For each hole compute its size.

    Holes should be colored as -1.
    """
    grid = np.copy(grid)
    hole_sizes = []

    hole_counter = 0
    hole_sizes = []
    while True:
        # find a point where grid[i,j] = -1 (a hole)
        ihs,jhs = np.where(grid == -1)
        if len(ihs) == 0:
            # no more holes
            break
        # pick the first hole
        i0,j0 = ihs[0], jhs[0]
        hole_counter += 1
        # fill thils hole with color 
        color = hole_counter + 1
        grid = grid_flood_fill(grid, i0, j0, target_color=-1, replacement_color=color)
        # count number of squares belonging to this hole
        ic,jc = np.where(grid == color)
        size = len(ic)
        hole_sizes.append(size)

    print Flake(grid)
        
    return hole_sizes
        
def detect_edges(grid):
    """
    detect the edges of a grid using Sobel operator

    Parameters
    ----------
    grid      :  2d numpy array, integers

    Returns
    -------
    edges     :  positions that are 0 in the grid but have neighbours != 0
    """
    # Find edges using Sobe operator
    nx,ny = grid.shape
    edges = np.zeros((nx,ny))
    for i in range(0, nx):
        for j in range(0, ny):
            # horizontal edges
            if i == 0:
                edges[i,j] +=               - 1*grid[i,j] + 1*grid[i+1,j]
            elif i == nx-1:
                edges[i,j] += 1*grid[i-1,j] - 1*grid[i,j]
            else:
                edges[i,j] += 1*grid[i-1,j] - 2*grid[i,j] + 1*grid[i+1,j]
            # vertical edges
            if j == 0:
                edges[i,j] +=               - 1*grid[i,j] + 1*grid[i,j+1]
            elif j == ny-1:
                edges[i,j] += 1*grid[i,j-1] - 1*grid[i,j]
            else:
                edges[i,j] += 1*grid[i,j-1] - 2*grid[i,j] + 1*grid[i,j+1]

    return edges
                
def test_random_selection():
    probs = np.array([0.2, 0.7, 0.1])
    counts = np.array([0.0, 0.0, 0.0])
    for it in range(0, 10000):
        i = random_selection(probs)
        counts[i] += 1

    observed_probs = counts/np.sum(counts)
    print "probabilities  %s" % probs
    print "observed probabilities  %s" % observed_probs
    assert la.norm(probs - observed_probs) < 0.05

def grow_flakes_monte_carlo(max_generations, connectivity="meso_beta"):
    """
    Grow a porphyrin flake for `max_generations` by attaching monomers at random 
    with equal probability at every posible position.
    """
    if connectivity in ["meso_meso_wire", "meso_beta_wire"]:
        size = (max_generations+1)*2
    else:
        size = max((max_generations+1)*2,11)
    gr = np.zeros((size,size), dtype=int)
    orient = np.zeros((size,size), dtype=int)
    # place a monomer in the center
    gr[size/2,size/2] = 1
    orient[size/2,size/2] = 0

    flake = Flake(gr, orient)
    for gen in range(1, max_generations+1):
        # List of candidate flakes that are bigger by one monomer
        next_flakes = flake.grow_connectivity(connectivity)
        # Randomly select one of them
        nr = len(next_flakes)
        # All candidates have the same probability 1/n
        probs = np.ones(nr)/float(nr)
        # index of selected flake
        i_sel = random_selection(probs)
        flake = next_flakes[i_sel]        
        # find holes, once produced they cannot be filled anymore
        flake.grid = identify_holes(flake.grid)

        #print flake

    print Flake(flake.grid, orient)

    # analyse holes in grown island
    sizes =  sizes_of_holes(flake.grid)
    print "number of holes : %d" % len(sizes)
    print "sizes of Holes : %s" % sizes
    
    return flake

def plot_polyomino_function(grid, values, ax):
    """
    Plot a function on a grid. values[i,j] are the function values at the grid
    position (i,j).
    """
    # colormap for translating values into colors
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap('cividis') #'seismic')

    ax.clear()
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect('equal', 'datalim')
    
    # plot a little square around each occupied grid point
    nx,ny = grid.shape
    for i in range(0, nx):
        for j in range(0, ny):
            if grid[i,j] == 1:
                # corners of square, starting in the lower left corner
                #  ___
                #  |x|
                #  ---
                # The last point has to be repeated to close the loop.
                xpos = np.array([i-0.5,i+0.5,i+0.5,i-0.5,i-0.5])
                ypos = np.array([j-0.5,j-0.5,j+0.5,j+0.5,j-0.5])
                # x- and y-axis are exchanged so that the flakes look the
                # same in the plots as on the terminal.
                ax.plot(ypos,-xpos, ls="-", lw=1.5, color="black", alpha=0.6)
                ax.fill(ypos,-xpos, color=cmap(values[i,j]))


def grow_flakes_monte_carlo_fast(max_generations,
                                 growth_probabilities="uniform",
                                 num_samples=10):
    """
    Grow a porphyrin flake for `max_generations` by attaching monomers at random 
    with equal probability at every posible position.

    Parameters
    ----------
    growth_probabilities  :    'uniform' - any position is equally likely to be filled
                               'orbital' - probability of attaching a new monomer at a 
                                           position depends on the a_2u orbital coefficients
                                           at neighbouring sites
    """
    ## matrix elements between Gouterman orbitals fused horizontally or vertically
    orbe, Sh,Sv, H0h,H0v = Gouterman_matelems()
    # The second Gouterman orbital has symmetry a2u, this is the orbital
    # that has electron density at the meso carbons.
    index_a2u = 1

    # Holds number of holes in each grown sample
    holes = []
    
    for s in range(0, num_samples):
        print "sample %4.1d/%4.1d"  % (s+1, num_samples)
    
        size = (max_generations+1)*2
        gr = np.zeros((size,size), dtype=int)
        
        gr = np.zeros((size,size), dtype=int)
        orient = np.zeros((size,size), dtype=int)
        # place a monomer in the center
        gr[size/2,size/2] = 1
        orient[size/2,size/2] = 0
    
        flake = Flake(gr)

        for gen in range(1, max_generations+1):
            if gen % 100 == 0:
                print "   generation %4.1d/%4.1d"  % (gen, max_generations)
            
            # Find edges where a new monomer can be attached
            grid_ = np.copy(flake.grid)
            # Set all positions which are not empty or outside the
            # flake to 1.
            grid_[grid_ != 0] = 1
            # 
            edges = detect_edges(grid_)
            # new flakes can be attached only at empty positions
            exterior_edges = np.copy(edges)
            exterior_edges[grid_ != 0] = 0
            
            print "Flake"
            print Flake(flake.grid)
            #print "Exterior Edges"
            #print Flake(exterior_edges)
            
            # growth_probs[i,j] is the probability that a new monomer
            # is attached at position (i,j)
            growth_probs = np.zeros(exterior_edges.shape)
            
            if growth_probabilities == "orbital":
                en_tot, HLgap, orbeHueckel, orbsHueckel = flake.hueckel(orbe, Sh,Sv, H0h,H0v)
                # a2u coefficients in HOMO
                orbH = orbsHueckel[:,flake.HOMO]
                coeffs_a2u = np.zeros(flake.grid.shape)
                # enumerate porphyrin site
                for k in range(0, flake.N):
                    i,j = flake.chain2grid[k]
                    # index of a2u orbital
                    idx = k * flake.norb + index_a2u
                    coeffs_a2u[i,j] = orbH[idx]
                    # sum of modulus squares of a2u coefficients
                prob_tot_a2u = np.sum(abs(coeffs_a2u)**2)
                print np.sum(abs(orbH)**2)
                #print "probability of a2u component in HOMO: %e" % prob_tot_a2u

                if gen > 50:
                    prob_a2u = abs(coeffs_a2u)**2
                    # plot sketch of polyomino and a2u coefficients
                    ax = plt.gca()
                    # normalize probability such that maximum equals 1
                    plot_polyomino_function(flake.grid, prob_a2u / prob_a2u.max(), ax)
                    plt.show()

                growth_probs[exterior_edges != 0] = prob_tot_a2u[exterior_edges != 0]
            
            else:
                assert growth_probabilities == "uniform"
                # all positions are equally possible
                growth_probs[exterior_edges != 0] = 1

            #print "Growth probabilities"
            #print Flake(growth_probs)
                
            # normalize probabilities
            growth_probs /= np.sum(growth_probs)
            
            probs = growth_probs.flatten()
            # index of selected growth position
            i_sel = random_selection(probs)
            # convert i_sel into grid position (i,j)
            nx, ny = growth_probs.shape
            i = i_sel / ny
            j = i_sel % ny

            # growth position should be empty
            assert flake.grid[i,j] == 0

            next_grid = np.copy(flake.grid)
            next_grid[i,j] = 1
            
            # next, bigger flake
            flake = Flake(next_grid, orient)
            
            # find holes, once produced they cannot be filled anymore
            flake.grid = identify_holes(flake.grid)

            
        # analyse holes in grown island
        sizes =  sizes_of_holes(flake.grid)
        print "number of holes : %d" % len(sizes)
        print "sizes of Holes : %s" % sizes

        num_holes = len(sizes)
        holes.append(num_holes)
        
    return holes


if __name__ == "__main__":
    #test_random_selection()
    
    max_generations = 100
    num_samples = 1000
    connectivity = 'meso_beta'
    
    #grow_flakes_monte_carlo(max_generations, connectivity)

    holes = grow_flakes_monte_carlo_fast(max_generations)

        
    
