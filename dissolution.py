import numpy as np
import scipy.sparse as spr

from config import simInputData
from utils import solve_equation


def create_vector(sid:simInputData, in_nodes):
    """ Creates vector result for B concentration calculation.
    
    For inlet nodes elements of the vector correspond explicitly
    to the concentration in nodes, for other it corresponds to 
    mixing condition.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        nsq - number of nodes in the network squared
        cb_in - B concentration in inlet nodes

    in_nodes : list
        list of inlet nodes   
    
    Returns
    ------
    scipy sparse vector
        result vector for B concentration calculation
    """
    data, row, col = [], [], []
    for node in in_nodes:
        data.append(sid.cb_in)
        row.append(node)
        col.append(0)
    return spr.csc_matrix((data, (row, col)), shape=(sid.nsq, 1))

def find_cb(sid, diams, lens, flow, inc_matrix, in_nodes, out_nodes, cb_b, triangles_inc, vol_a):
    if sid.include_vol_a:
        return find_cb_vol(sid, diams, lens, flow, inc_matrix, in_nodes, out_nodes, cb_b, triangles_inc, vol_a)
    else:
        return find_cb_novol(sid, diams, lens, flow, inc_matrix, in_nodes, out_nodes, cb_b)


def find_cb_novol(sid, diams, lens, flow, inc_matrix, in_nodes, out_nodes, cb_b):
    """ Calculates B concentration.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        Da - Damkohler number
        G - diffusion-reaction relation coefficient

    diams : numpy array
        current diameters of edges
    lens : numpy array
        lengths of edges
    flow : numpy array
        vector of flows
    inc_matrix : scipy sparse array
        incidence matrix
    in_nodes : list
        list of inlet nodes
    out_nodes : list
        list of outlet nodes
    cb_b : scipy sparse vector
        result vector for concentration calculation
    
    Returns
    -------
    cb : numpy array
        vector of B concentration in nodes
    """
    cb_inc = np.abs(inc_matrix.transpose() @ (spr.diags(flow) @ inc_matrix > 0)) # find incidence for cb (only upstream flow matters)
    qc = flow * np.exp(-sid.Da / (1 + sid.G * diams) * diams * lens / np.abs(flow)) # find vector with non-diagonal coefficients
    qc_matrix = np.abs(inc_matrix.transpose() @ spr.diags(qc) @ inc_matrix)
    cb_matrix = cb_inc.multiply(qc_matrix)
    diag = -np.abs(inc_matrix.transpose()) @ np.abs(flow) / 2 # find diagonal coefficients (inlet flow for each node)
    for node in in_nodes:
        diag[node] = 1 # set diagonal for input nodes to 1
    for node in out_nodes:
        if diag[node] != 0: # fix for nodes which are connected only to other out_nodes - without it we get a singular matrix (whole row of zeros)
            diag[node] *= 2 # multiply diagonal for output nodes (they have no outlet, so inlet flow is equal to whole flow)
        else:
            diag[node] = 1
    cb_matrix.setdiag(diag) # replace diagonal
    cb = solve_equation(cb_matrix, cb_b)
    return cb

def find_cb_vol(sid, diams, lens, flow, inc_matrix, in_nodes, out_nodes, cb_b, triangles_inc, vol_a):
    """ Calculates B concentration.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        Da - Damkohler number
        G - diffusion-reaction relation coefficient

    diams : numpy array
        current diameters of edges
    lens : numpy array
        lengths of edges
    flow : numpy array
        vector of flows
    inc_matrix : scipy sparse array
        incidence matrix
    in_nodes : list
        list of inlet nodes
    out_nodes : list
        list of outlet nodes
    cb_b : scipy sparse vector
        result vector for concentration calculation
    
    Returns
    -------
    cb : numpy array
        vector of B concentration in nodes
    """
    cb_inc = np.abs(inc_matrix.transpose() @ (spr.diags(flow) @ inc_matrix > 0)) # find incidence for cb (only upstream flow matters)
    vol_a_avail = (triangles_inc @ vol_a > 0) * 1
    #qc = flow * np.exp(-sid.Da / (1 + sid.G * diams) * diams * lens / np.abs(flow)) # find vector with non-diagonal coefficients
    diag = -np.abs(inc_matrix.transpose()) @ np.abs(flow) / 2 # find diagonal coefficients (inlet flow for each node)
    for node in in_nodes:
        diag[node] = 1 # set diagonal for input nodes to 1
    for node in out_nodes:
        if diag[node] != 0: # fix for nodes which are connected only to other out_nodes - without it we get a singular matrix (whole row of zeros)
            diag[node] *= 2 # multiply diagonal for output nodes (they have no outlet, so inlet flow is equal to whole flow)
        else:
            diag[node] = 1    

    exp_eff0 = -1
    exp_eff  = np.exp(-sid.Da / (1 + sid.G * diams) * diams * lens / np.abs(flow))
    exp_eff00 = exp_eff.copy()
    edge_triangles = np.array(np.sum(triangles_inc, axis = 1))[:, 0] # calculate how many triangles each edge has as neighbors (1 or 2)
    it_k = 0
    alfa = 1
    while np.linalg.norm(exp_eff - exp_eff0) > sid.k_it_th:
        print (f'K iterations {it_k}')
        print (np.linalg.norm(exp_eff - exp_eff0))
        it_k += 1
        exp_eff0 = exp_eff.copy()
        qc = flow * vol_a_avail * exp_eff + flow * (1 - vol_a_avail)
        # if there is available volume, we include reduction of cb in pore, if not, then q_in cb_in = q_out cb_out
        qc_matrix = np.abs(inc_matrix.transpose() @ spr.diags(qc) @ inc_matrix)
        cb_matrix = cb_inc.multiply(qc_matrix)
        cb_matrix.setdiag(diag) # replace diagonal
        cb = solve_equation(cb_matrix, cb_b)

        growth_matrix = np.abs((spr.diags(flow) @ inc_matrix > 0))
        cb_growth = growth_matrix @ cb
        growth = cb_growth * np.abs(flow) * 2 / sid.Da * (1 - exp_eff00) * vol_a_avail * sid.dt
        #vol_a_dissolved = triangles_inc.T @ (cb_growth * np.abs(flow) * 2 / sid.Da * (1 - exp_eff) * vol_a_avail / 2 * sid.dt)
        vol_a_dissolved = triangles_inc.T @ growth / 2
        vol_rat = np.minimum(1, vol_a / vol_a_dissolved)
        vol_rat[np.isnan(vol_rat)] = 1
        alfa = triangles_inc @ vol_rat / edge_triangles
        exp_eff = exp_eff00 * alfa
    print (np.linalg.norm(exp_eff - exp_eff0))
    return cb, alfa
