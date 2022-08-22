from collections import defaultdict
import numpy as np
import scipy.sparse as spr

from utils import d_update

from config import simInputData


def create_vector(sid:simInputData, in_nodes, out_nodes):
    presult = np.zeros(sid.nsq)
    for node in in_nodes:
        presult[node] = sid.pin
    for node in out_nodes:
        presult[node] = sid.pout
    return presult

def update_matrix(sid:simInputData, edges, in_nodes, out_nodes):

    data, row, col = [], [], []
    insert = defaultdict(float)
    diag = np.zeros(sid.nsq)
    for n1, n2, d, l, t in edges:
        res = sid.c1 / sid.mu * d ** 4 / l
        if t == 0:
            data.append(res)
            row.append(n1)
            col.append(n2)
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n1] -= res
            diag[n2] -= res
        elif t == 1:
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
        elif t == 2:
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
    for node, datum in enumerate(diag):
        if datum != 0:
            row.append(node)
            col.append(node)
            data.append(datum)

    for node in in_nodes+out_nodes:
        row.append(node)
        col.append(node)
        data.append(1)

    return spr.csr_matrix((data, (row, col)), shape=(sid.nsq, sid.nsq))

def update_graph(sid:simInputData, edges, pnow):
    d_pres = np.zeros(len(edges))
    if sid.shear_d:
        for i,e in enumerate(edges):
            n1, n2, d, l, t = e
            F = d * np.abs(pnow[n1] - pnow[n2]) / 2 / l
            d_pres[i] = d_update(F, sid.F_p)
    return d_pres


def update_network(G1, sid:simInputData, edges, pnow):
    Q_in = 0
    Q_out = 0
    in_edges = 0

    for n1, n2, d, l, t in edges:
        G1[n1][n2]['d']= d
        q = sid.c1 / sid.mu * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
        G1[n1][n2]['q'] = q
        
        if t == 1:    
            Q_in += q
            in_edges += 1
        if t == 2:
            Q_out += q
    
    print('Q_in =', Q_in, 'Q_out =', Q_out)

    G = sid.k * sid.noise[1] / sid.D / sid.alpha
    q_in = Q_in / in_edges
    Da = np.pi * sid.noise[1] * sid.k * sid.l / q_in / (1 + G)
    print ('G = ', G, ', Da = ', Da)
    return G1
