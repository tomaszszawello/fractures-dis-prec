import numpy as np
import scipy.sparse as spr

from config import simInputData

def create_matrices(sid:simInputData, edges, in_nodes, out_nodes):
    """ Creates incidence matrices for flow and concentration calculation.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        nsq - number of nodes in the network squared

    edges : list
        list of all edges in the network, containing:
        n1 - start node
        n2 - end node
        d - diameter
        l - length
        t - type (0 - regular, 1 - regular to input, 2 - regular to output)

    in_nodes : list
        list of input nodes
    
    out_nodes : list
        list of output nodes

    Returns
    -------
    sparse csr matrix
        matrix of incidence (ne x nsq)
    
    sparse csr matrix
        matrix of middle nodes (nsq x nsq)
    
    sparse csr matrix
        diagonal matrix of input and output nodes (nsq x nsq)
    
    sparse csr matrix
        matrix of connections for inlet edges (ne x nsq)

    diams : list
        vector of edge diameters
    
    lens : list
        vector of edge lengths

    """
    ne = len(edges)
    data, row, col = [], [], [] # data for standard incidence matrix (ne x nsq)
    diams, lens = [], [] # vectors of diameters and lengths
    data_mid, row_mid, col_mid = [], [], [] # data for matrix keeping connections of only middle nodes (nsq x nsq)
    data_bound, row_bound, col_bound = [], [], [] # data for diagonal matrix for input and output (nsq x nsq)
    data_in, row_in, col_in = [], [], [] # data for matrix keeping connections of only input nodes
    reg_nodes = [] # list of regular nodes
    for i, e in enumerate(edges):
        n1, n2, d, l, t = e
        data.append(-1)
        row.append(i)
        col.append(n1)
        data.append(1)
        row.append(i)
        col.append(n2)
        diams.append(d)
        lens.append(l)
        if t == 0:
            data_mid.extend((1, 1))
            row_mid.extend((n1, n2))
            col_mid.extend((n2, n1))
            reg_nodes.extend((n1, n2))
        elif t == 1:
            data_mid.append(1)
            row_mid.append(n1)
            col_mid.append(n2)
            data_in.append(1)
            row_in.append(i)
            col_in.append(n1)
            data_in.append(-1)
            row_in.append(i)
            col_in.append(n2)
    for node in in_nodes + out_nodes:
        data_bound.append(1)
        row_bound.append(node)
        col_bound.append(node)
    reg_nodes = list(set(reg_nodes))
    for node in reg_nodes:
        data_mid.append(1)
        row_mid.append(node)
        col_mid.append(node)
    return spr.csr_matrix((data, (row, col)), shape=(ne, sid.nsq)), \
        spr.csr_matrix((data_mid, (row_mid, col_mid)),shape=(sid.nsq, sid.nsq)), \
        spr.csr_matrix((data_bound, (row_bound, col_bound)), shape=(sid.nsq, sid.nsq)), \
        spr.csr_matrix((data_in, (row_in, col_in)), shape=(ne, sid.nsq)), \
        np.array(diams), np.array(lens)
