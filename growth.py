import numpy as np
import scipy.sparse as spr

from config import simInputData

def update_apertures(sid, flow, cb, apertures, lens, inc_matrix, out_edges, dt):
    """ Updates diameters in case of dissolution.

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

    out_edges : numpy array
        vector with 1 for outlet edges and 0 otherwise

    Returns
    -------
    diams : numpy array
        new diameters of edges

    dt : float
        current timestep
    
    breakthrough : bool
        information whether network is dissolved
    """
    breakthrough = False
    cb_growth = np.abs((spr.diags(flow) @ inc_matrix > 0)) @ cb # creates list of concentrations which should be used for growth of each edge (upstream one)
    aperture_change = cb_growth * np.abs(flow)  / (sid.Da * lens) * (1 - np.exp(-2 * sid.Da / (1 + sid.G * apertures) * lens / np.abs(flow)))
    if sid.adaptive_dt:
        dt = sid.growth_rate / np.max(aperture_change / apertures)
        if dt > sid.dt_max:
            dt = sid.dt_max
    else:
        dt = sid.dt
    apertures += dt * aperture_change
    if np.max(out_edges * apertures) > sid.b_break:
        breakthrough = True
        print ('Network dissolved.')
    return apertures, dt, breakthrough