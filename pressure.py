from collections import defaultdict
import numpy as np
import networkx as nx
import scipy.sparse as spr

from utils import d_update

from config import simInputData


def create_vector(sid:simInputData, in_nodes, out_nodes, edges):
    # in_edges = 0
    # for n1, n2, d, l, t in edges:
    #     if t == 1:\\\\\\\\\\\\\\\\\
    #         in_edges += 1
    data, row, col = [], [], []
    presult = np.zeros(sid.nsq)
    for node in in_nodes:
        data.append(1)
        row.append(node)
        col.append(0)
    #for node in out_nodes:
    #    data.append(-sid.qin * 2 * len(in_nodes))
    #    row.append(node)
    #    col.append(0)
    return spr.csc_matrix((data, (row, col)), shape=(sid.nsq, 1))



def update_matrix(sid:simInputData, edges, in_nodes, out_nodes):

    data, row, col = [], [], []
    insert = defaultdict(float)
    diag = np.zeros(sid.nsq)
    for n1, n2, d, l, t in edges:
        res = d ** 4 / l
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
            insert[n1] += res
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

    for node in out_nodes:
        row.append(node)
        col.append(node)
        data.append(1)

    sum_insert = sum(insert.values())

    for node in in_nodes:
        data.append(-sum_insert)
        row.append(node)
        col.append(node)

        for ins_node, ins in insert.items():
            data.append(ins)
            row.append(node)
            col.append(ins_node)

    return spr.csr_matrix((data, (row, col)), shape=(sid.nsq, sid.nsq))



def update_network(G1, sid:simInputData, edges, diams, flows, in_nodes, out_nodes):
    Q_in = 0
    Q_out = 0
    
    #nx.set_edge_attributes(G1, diams, name="d")

#    for i, e in enumerate(G1.edges()):
#        n1, n2 = e
    for i, e in enumerate(edges):
        n1, n2, d, l, t = e
        d = diams[i]
        q = flows[i]
        G1[n1][n2]['d']= d
        G1[n1][n2]['q'] = q
        
        if n1 in in_nodes or n2 in in_nodes:    
            Q_in += np.abs(q)
        if n1 in out_nodes or n2 in out_nodes:
            Q_out += np.abs(q)

    print('Q_in =', Q_in, 'Q_out =', Q_out)

    return G1
