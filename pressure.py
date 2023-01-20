import networkx as nx
import numpy as np
import scipy.sparse as spr

from config import simInputData
from utils import solve_equation


def create_vector(sid:simInputData, in_nodes):
    """ Creates vector result for pressure calculation.
    
    For inlet and outlet nodes elements of the vector correspond explicitly
    to the pressure in nodes, for regular nodes elements of the vector equal
    0 correspond to flow continuity.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        nsq - number of nodes in the network squared

    in_nodes : list
        list of inlet nodes
    
    Returns
    -------
    scipy sparse vector
        result vector for pressure calculation
    """
    data, row, col = [], [], []
    for node in in_nodes:
        data.append(1)
        row.append(node)
        col.append(0)
    return spr.csc_matrix((data, (row, col)), shape=(sid.nsq, 1))

def find_flow(sid, apertures, lens, inc_matrix, mid_matrix, bound_matrix, in_matrix, pressure_b, in_edges):
    """ Calculates pressure and flow.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        qin - characteristic flow for inlet edge

    diams : numpy array
        current diameters of edges

    lens : numpy array
        lengths of edges

    inc_matrix : scipy sparse array
        incidence matrix

    mid_matrix : scipy sparse array
        matrix zeroing rows for input and output nodes

    bound_matrix : scipy sparse array
        diagonal matrix with ones for input and output nodes

    in_matrix : scipy sparse array
        incidence matrix for inlet edges

    pressure_b : scipy sparse vector
        result vector for pressure equation

    in_nodes : list
        list of inlet nodes


    Returns
    -------
    pressure : numpy array
        vector of pressure in nodes

    flow : numpy array
        vector of flows in edges
    """
    p_matrix = inc_matrix.transpose() @ spr.diags(apertures ** 3 / lens) @ inc_matrix # create matrix
    p_matrix = p_matrix.multiply(mid_matrix) + bound_matrix
    pressure = solve_equation(p_matrix, pressure_b)
    q_in = np.abs(np.sum(apertures ** 3 / lens * (in_matrix @ pressure))) # calculate inlet flow
    pressure *= sid.qin * np.sum(in_edges) / q_in # normalize pressure to match condition for constant inlet flow
    flow = apertures ** 3 / lens * (inc_matrix @ pressure)
    return pressure, flow


def update_network(G, edge_list, apertures, flow):
    nx.set_edge_attributes(G, dict(zip(edge_list, apertures)), 'b')
    nx.set_edge_attributes(G, dict(zip(edge_list, flow)), 'q')
    return G