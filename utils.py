import numpy as np
import os
import scipy.sparse.linalg as sprlin
import scipy.sparse as spr

def solve_equation(A, b):
    """ Solves matrix equation A * x = b.

    Parameters
    -------
    A : scipy sparse matrix
        matrix A from equation

    b : scipy sparse vector
        result b from equation

    Returns
    -------
    scipy sparse vector
        result x from equation    
    """
    return sprlin.spsolve(A, b)

def initialize_iterators(sid):
    """ Creates iterators for simulation steps, time and other conditions.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        iters - max no. of iterations of new simulation
        tmax - max time of new simulation
        old_iters - no. of previous iterations (if loaded from saved file)
        old_t - time of previous simulation (if loaded from saved file)
        dt - initial timestep

    Returns
    -------
    iters : int 
        max no. of new iterations

    i : int
        iterator in range from old iterations to sum of old and new

    tmax : float
        max new time

    t : float
        time iterator in range from old time to sum of old and new

    dt : float
        initial timestep

    breakthrough : bool
        parameter stating if the system was dissolved (if diameter of edge
        going to the output grew at least sid.d_break times)
    """
    iters = sid.old_iters + sid.iters
    tmax = sid.old_t + sid.tmax
    i = sid.old_iters
    t = sid.old_t
    dt = sid.dt
    breakthrough = False
    return iters, tmax, i, t, dt, breakthrough

def update_iterators(sid, i, t, dt):
    """ Updates iterators in simulation and in configuration class.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        old_iters - no. of previous iterations (if loaded from saved file)
        old_t - time of previous simulation (if loaded from saved file)

    i : int
        current iteration

    t : float
        current time

    dt : float
        current timestep

    Returns
    -------
    i : int
        current iteration

    t : float
        current time
    """
    i += 1
    sid.old_iters += 1 # update simulation iterations in configuration class 
    t += dt
    sid.old_t += dt # update simulation time in configuration class
    return i, t



def make_dir(sid):
        i = 0
        dirname2 = sid.dirname
        while (sid.dirname == dirname2):
            if not os.path.isdir(sid.dirname + "/" + str(i)):
                sid.dirname = sid.dirname + "/" + str(i)
            else:
                i += 1
        os.makedirs(sid.dirname)