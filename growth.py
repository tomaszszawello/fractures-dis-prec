import numpy as np
import scipy.sparse as spr

from config import simInputData

def update_diameters(sid, flow, cb, cc, diams, lens, inc_matrix, triangles_inc, vols, out_edges, alfa):
    if sid.include_cc:
        if sid.include_vol_a:
            return update_diameters_dpv(sid, flow, cb, cc, diams, lens, inc_matrix, triangles_inc, vols, out_edges, alfa)
        else:
            return update_diameters_dp(sid, flow, cb, cc, diams, lens, inc_matrix, out_edges)
    else:
        return update_diameters_d(sid, flow, cb, diams, lens, inc_matrix, out_edges)
        

def update_diameters_d(sid, flow, cb, diams, lens, inc_matrix, out_edges):
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
    cb_growth = np.abs((spr.diags(flow) @ inc_matrix > 0)) @ cb # creates list of concentrations which should be used for growth of each edge (upstream one)
    diameter_change = cb_growth * np.abs(flow)  / (sid.Da * lens * diams) * (1 - np.exp(-sid.Da / (1 + sid.G * diams) * diams * lens / np.abs(flow)))
    if sid.adaptive_dt:
        dt = sid.growth_rate / np.max(diameter_change / diams)
        if dt > sid.dt_max:
            dt = sid.dt_max
    else:
        dt = sid.dt
    diams += dt * diameter_change
    if np.max(out_edges * diams) > sid.d_break:
        breakthrough = True
        print ('Network dissolved.')
    return diams, None, dt, breakthrough

def update_diameters_dp(sid, flow, cb, cc, diams, lens, inc_matrix, out_edges):
    """ Updates diameters in case of dissolution + precipitation.

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
    
    breakthrough : bool
        information whether network is dissolved
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
    diams = (diams > sid.dmin) * diams + (diams <= sid.dmin) * sid.dmin
    if np.max(out_edges * diams) > sid.d_break:
        breakthrough = True
        print ('Network dissolved.')
    return diams, None, dt, breakthrough

def update_diameters_dpv(sid, flow, cb, cc, diams, lens, inc_matrix, triangles_inc, vols, out_edges, alfa):
    """ Updates diameters in case of dissolution + precipitation.

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
    
    breakthrough : bool
        information whether network is dissolved
    """
    # create list of concentrations which should be used for
    # growth of each edge (upstream one)
    breakthrough = False
    growth_matrix = np.abs((spr.diags(flow) @ inc_matrix > 0))
    cb_growth = growth_matrix @ cb
    cc_growth = growth_matrix @ cc
    vol_a_avail = (triangles_inc @ vols > 0) * 1
    diameter_growth = alfa * cb_growth * np.abs(flow)  / (sid.Da * lens * diams) * (1
    - np.exp(-sid.Da / (1 + sid.G * diams) * diams * lens / np.abs(flow))) * vol_a_avail
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
    #vols -= triangles_inc.T @ (diams * diameter_growth * dt * lens / 2)
    vols -= triangles_inc.T @ (alfa * cb_growth * np.abs(flow) * 2 / sid.Da * (1 - np.exp(-sid.Da / (1 + sid.G * diams) * diams * lens / np.abs(flow))) * vol_a_avail / 2 * sid.dt)
    #vols -= triangles_inc.T @ diameter_growth * sid.dt / 2
    diams += dt * diameter_change
    diams = (diams > sid.dmin) * diams + (diams <= sid.dmin) * sid.dmin
    #vols -= triangles_inc.T @ (diameter_growth * dt / 2)
    vols = (vols > 0) * vols
    if np.max(out_edges * diams) > sid.d_break:
        breakthrough = True
        print ('Network dissolved.')
    return diams, vols, dt, breakthrough