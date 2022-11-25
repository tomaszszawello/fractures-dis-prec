import numpy as np
import scipy.sparse as spr

from config import simInputData
from utils import solve_equation


def create_vector(sid, diams, lens, flow, inc_matrix, in_nodes, cb):
    """ Creates vector result for C concentration calculation.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        Da - Damkohler number
        G - diffusion-reaction relation coefficient
        K - ratio of precipitation rate to dissolution rate
        cc_in - C concentration in inlet nodes

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

    cb : scipy sparse vector
        vector of B concentration
    

    Returns
    -------
    cc : numpy array
        vector result for C concentration calculation
    """
    # find incidence for cb (only upstream flow matters)
    cb_inc = np.abs(inc_matrix.transpose() @ (spr.diags(flow) @ inc_matrix > 0))
    # find vector with non-diagonal coefficients
    qc = flow / (sid.K - 1) * (np.exp(-sid.Da / (1 + sid.G * diams) * diams
    * lens / np.abs(flow)) - np.exp(-sid.Da * sid.K / (1 + sid.G * sid.K
    * diams) * diams * lens / np.abs(flow)))
    qc_matrix = np.abs(inc_matrix.transpose() @ spr.diags(qc) @ inc_matrix)
    cb_matrix = cb_inc.multiply(qc_matrix)
    cb_matrix.setdiag(np.zeros(sid.nsq)) # set diagonal to zero
    cc_b = -cb_matrix @ cb
    for node in in_nodes:
        cc_b[node] = sid.cc_in # set result for input nodes to cc_in
    return cc_b


def find_cc(sid, diams, lens, flow, inc_matrix, in_nodes, out_nodes, cc_b):
    """ Calculates pressure and flow.

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

    cc_b : scipy sparse vector
        result vector for concentration calculation

    cb : numpy array
        vector of B concentration in nodes
    

    Returns
    -------
    cc : numpy array
        vector of C concentration in nodes
    """
    # find incidence for cc (only upstream flow matters)
    cc_inc = np.abs(inc_matrix.transpose() @ (spr.diags(flow) @ inc_matrix > 0))
    # find vector with non-diagonal coefficients
    qc = flow * np.exp(-sid.Da * sid.K / (1 + sid.G * sid.K * diams) * diams
    * lens / np.abs(flow))
    qc_matrix = np.abs(inc_matrix.transpose() @ spr.diags(qc) @ inc_matrix)
    cc_matrix = cc_inc.multiply(qc_matrix)
    # find diagonal coefficients (inlet flow for each node)
    diag = -np.abs(inc_matrix.transpose()) @ np.abs(flow) / 2
    for node in in_nodes:
        diag[node] = 1 # set diagonal for input nodes to 1
    for node in out_nodes:
        if diag[node] != 0: # fix for nodes which are connected only to other 
            # out_nodes - without it we get a singular matrix (whole row of zeros)
            diag[node] *= 2 # multiply diagonal for output nodes
            # (they have no outlet, so inlet flow is equal to whole flow)
        else:
            diag[node] = 1
        diag[node] *= 2 
    cc_matrix.setdiag(diag) # replace diagonal
    cc = solve_equation(cc_matrix, cc_b)
    return cc
