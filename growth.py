""" Updates edges diameters based on dissolution and precipitation.

This module calculates the change od diameters in the network, resulting from
dissolution (and precipitation, if enabled). Based on that change, new
timestep is calculated.

Notable functions
-------
update_diameters(SimInputData, Incidence, Edges, np.ndarray, np.ndarray) \
    -> tuple[bool, float]
    update diameters, calculate timestep and check if network is dissolved
"""

import numpy as np
import scipy.sparse as spr

from config import SimInputData
from incidence import Edges, Incidence


def update_apertures(sid: SimInputData, inc: Incidence, edges: Edges, \
    cb: np.ndarray) -> tuple[bool, float]:
    """ Update diameters.

    This function updates diameters of edges, calculates the next timestep (if
    adt is used) and checks if the network is dissolved. Based on config, we
    include either dissolution or both dissolution and precipitation.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation
        include_cc : bool
        dmin : float
        dmin_th : float
        d_break : float
        include_adt : bool
        growth_rate : float
        dt : float
        dt_max : float

    inc : Incidence class object
        matrices of incidence

    edges : Edges class object
        all edges in network and their parameters
        diams : numpy ndarray (ne)
        lens : numpy ndarray (ne)
        outlet : numpy ndarray (ne)

    cb : numpy ndarray (nsq)
        vector of substance B concentration

    cc : numpy ndarray (nsq)
        vector of substance C concentration

    Returns
    -------
    breakthrough : bool
        parameter stating if the system was dissolved (if diameter of output
        edge grew at least to sid.d_break)

    dt_next : float
        new timestep
    """

    breakthrough = False
    # create list of concentrations which should be used for growth of each
    # edge (upstream one)
    cb_in = np.abs((spr.diags(edges.flow) @ inc.incidence > 0)) @ cb 
    aperture_change = cb_in * np.abs(edges.flow)  / (sid.Da * edges.lens) * (1 \
        - np.exp(-np.abs(2 * sid.Da / (1 + sid.G * edges.apertures) \
        * edges.lens / edges.flow)))
    if sid.include_adt:
        dt = sid.growth_rate / np.max(aperture_change / edges.apertures)
        if dt > sid.dt_max:
            dt = sid.dt_max
    else:
        dt = sid.dt
    edges.apertures += dt * aperture_change
    if np.max(edges.outlet * edges.apertures) > sid.b_break and sid.include_breakthrough:
        breakthrough = True
        print ('Network dissolved.')
    return breakthrough, dt
