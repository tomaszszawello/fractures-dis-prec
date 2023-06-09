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

def find_cb(sid, apertures, fracture_lens, lens, flow, inc_matrix, in_nodes, out_nodes, cb_b):
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
    qc = fracture_lens * flow * np.exp(-2 * sid.Da / (1 + sid.G * apertures) * lens / np.abs(flow)) # find vector with non-diagonal coefficients
    qc_matrix = np.abs(inc_matrix.transpose() @ spr.diags(qc) @ inc_matrix)
    cb_matrix = cb_inc.multiply(qc_matrix)
    diag = -np.abs(inc_matrix.transpose()) @ np.abs(flow * fracture_lens) / 2 # find diagonal coefficients (inlet flow for each node)
    for node in in_nodes:
        diag[node] = 1 # set diagonal for input nodes to 1
    for node in out_nodes:
        diag[node] *= 2 # multiply diagonal for output nodes (they have no outlet, so inlet flow is equal to whole flow)
    diag = diag * (diag != 0) + 1 * (diag == 0) # fix for nodes with 0 flow - without it we get a singular matrix
    cb_matrix.setdiag(diag) # replace diagonal
    cb = solve_equation(cb_matrix, cb_b)
    return cb