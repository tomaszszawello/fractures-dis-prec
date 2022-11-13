from smtplib import OLDSTYLE_AUTH
import networkx as nx
import numpy as np
import os
import scipy.sparse.linalg as sprlin
import scipy.sparse as spr
from delaunay import find_node

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

    breakthrough : bool
        parameter stating if the system was dissolved (if diameter of edge
        going to the output grew at least sid.d_break times)
    """
    iters = sid.old_iters + sid.iters
    tmax = sid.old_t + sid.tmax
    i = sid.old_iters
    t = sid.old_t
    breakthrough = False
    return iters, tmax, i, t, breakthrough

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

def update_diameters(sid, flow, cb, diams, lens, inc_matrix):
    """ Updates diameters.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        Da - Damkohler number
        G - 
        adaptive_dt - use adaptive timestep
        growth_rate - if adaptive timestep is on, maximum percentage of diameter which edge can grow
        dt_max - if adaptive timestep is on, maximal timestep
        dt - timestep if adaptive timestep is off

    flow : numpy array
        current flow in edges

    cb : numpy array
        current concentration in nodes

    diams : numpy array
        current diameters of edges

    lens : numpy array
        lengths of edges

    inc_matrix : scipy sparse array
        matrix of connections

    Returns
    -------
    diams : numpy array
        new diameters of edges

    dt : float
        current timestep
    """
    cb_growth = np.abs((spr.diags(flow) @ inc_matrix > 0)) @ cb # creates list of concentrations which should be used for growth of each edge (upstream one)
    diameter_change = cb_growth * np.abs(flow)  / (sid.Da * lens * diams) * (1 - np.exp(-sid.Da / (1 + sid.G * diams) * diams * lens / np.abs(flow)))
    if sid.adaptive_dt:
        dt = sid.growth_rate / np.max(diameter_change / diams)
        if dt > sid.dt_max:
            dt = sid.dt_max
    else:
        dt = sid.dt
    diams += dt * diameter_change
    return diams, dt

def update_diameters_pi(sid, flow, cb, cc, diams, lens, inc_matrix, out_edges):
    """ Updates diameters.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        Da - Damkohler number
        G - 
        adaptive_dt - use adaptive timestep
        growth_rate - if adaptive timestep is on, maximum percentage of diameter which edge can grow
        dt_max - if adaptive timestep is on, maximal timestep
        dt - timestep if adaptive timestep is off

    flow : numpy array
        current flow in edges

    cb : numpy array
        current B concentration in nodes

    cc : numpy array
        current C concentration in nodes

    diams : numpy array
        current diameters of edges

    lens : numpy array
        lengths of edges

    inc_matrix : scipy sparse array
        matrix of connections

    out_edges : numpy array
        vector with 1 for outlet edges and 0 otherwise

    Returns
    -------
    diams : numpy array
        new diameters of edges

    dt : float
        current timestep
    """
    # create list of concentrations which should be used for
    # growth of each edge (upstream one)
    breakthrough = False
    growth_matrix = np.abs((spr.diags(flow) @ inc_matrix > 0))
    cb_growth = growth_matrix @ cb
    cc_growth = growth_matrix @ cc
    diameter_growth = cb_growth * np.abs(flow)  / (sid.Da * lens * diams) * (1 
    - np.exp(-sid.Da / (1 + sid.G * diams) * diams * lens / np.abs(flow)))
    diameter_shrink_cb = cb_growth * np.abs(flow)  / (sid.Da * lens * diams
    * sid.Gamma) / (sid.K - 1) * (sid.K * (1 - np.exp(-sid.Da / (1 + sid.G
    * diams) * diams * lens / np.abs(flow))) - (1 - np.exp(-sid.Da * sid.K / (1
    + sid.G * sid.K * diams) * diams * lens / np.abs(flow))))
    diameter_shrink_cc = cc_growth * np.abs(flow)  / (sid.Da * lens * diams
    * sid.Gamma) * (1 - np.exp(-sid.Da * sid.K / (1 + sid.G * sid.K * diams)
    * diams * lens / np.abs(flow)))
    diameter_change = diameter_growth - diameter_shrink_cb - diameter_shrink_cc
    if sid.adaptive_dt:
        dt = sid.growth_rate / np.max(np.abs(diameter_change / diams))
        if dt > sid.dt_max:
            dt = sid.dt_max
    else:
        dt = sid.dt
    diams += dt * diameter_change
    diams = (diams > 0) * diams
    if np.max(out_edges * diams) > sid.d_break:
        breakthrough = True
    return diams, dt, breakthrough

def check_flow(flow, in_edges, out_edges):
    Q_in = np.sum(in_edges * np.abs(flow))
    Q_out = np.sum(out_edges * np.abs(flow))
    print('Q_in =', Q_in, 'Q_out =', Q_out)

def make_dir(sid):
        i = 0
        dirname2 = sid.dirname
        while (sid.dirname == dirname2):
            if not os.path.isdir(sid.dirname + "/" + str(i)):
                sid.dirname = sid.dirname + "/" + str(i)
            else:
                i += 1
        os.makedirs(sid.dirname)