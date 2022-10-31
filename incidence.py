import numpy as np
import scipy.sparse as spr

from config import simInputData

def create_incidence_matrix(sid:simInputData, edges):
    data, row, col = [], [], []
    data_d, row_d, col_d = [], [], []
    data_l, row_l, col_l = [], [], []
    diams = []
    lens = []
    ne = len(edges)
    for i, e in enumerate(edges):
        n1, n2, d, l, t = e
        if t == 0:
            data.append(1)
            row.append(i)
            col.append(n1)
            data.append(-1)
            row.append(i)
            col.append(n2)
        elif t == 1:
            data.append(1)
            row.append(i)
            col.append(n1)
            data.append(-1)
            row.append(i)
            col.append(n2)
        elif t == 2:
            data.append(1)
            row.append(i)
            col.append(n2)
            data.append(-1)
            row.append(i)
            col.append(n1)
        data_d.append(d)
        row_d.append(i)
        col_d.append(i)
        data_l.append(l)
        row_l.append(i)
        col_l.append(i)
        
    return spr.csc_matrix((data, (row, col)), shape=(ne, sid.nsq)), spr.csr_matrix((data_d, (row_d, col_d)), shape=(ne, ne)), spr.linalg.inv(spr.csr_matrix((data_l, (row_l, col_l)), shape=(ne, ne)))

def create_incidence_pressure_matrix(sid:simInputData, edges):
    data, row, col = [], [], []
    ne = len(edges)
    in_edges = []
    in_nodes = []
    out_edges = []
    out_nodes = []
    for i, e in enumerate(edges):
        n1, n2, d, l, t = e
        if t == 0:
            data.append(1)
            row.append(i)
            col.append(n1)
            data.append(-1)
            row.append(i)
            col.append(n2)
        elif t == 1:
            in_edges.append(i)
            in_nodes.append(n2)
            data.append(1)
            row.append(i)
            col.append(n1)
        elif t == 2:
            out_edges.append(i)
            out_nodes.append(n2)
            data.append(-1)
            row.append(i)
            col.append(n1)
    in_nodes = list(set(list(in_nodes)))
    out_nodes = list(set(list(out_nodes)))
    for edge in in_edges:
        for node in in_nodes:
                data.append(1)
                row.append(edge)
                col.append(node)
    for edge in out_edges:
        for node in out_nodes:
            data.append(1)
            row.append(edge)
            col.append(node)

    return spr.csr_matrix((data, (row, col)), shape=(ne, sid.nsq))

    
def create_in_out_matrix(sid:simInputData, edges, in_nodes, out_nodes):
    data, row, col = [], [], []
    ne = len(edges)
    in_edges = []
    in_nodes2 = []
    out_edges = []
    nodes = []
    data_vec_in, row_vec_in, col_vec_in = [], [], []
    data_in, row_in, col_in = [], [], []
    data_diag, row_diag, col_diag = [], [], []
    for i, e in enumerate(edges):
        n1, n2, d, l, t = e
        if t == 0:
            data.append(1)
            row.append(i)
            col.append(n1)
            data.append(-1)
            row.append(i)
            col.append(n2)
            data_in.append(1)
            row_in.append(n1)
            col_in.append(n2)
            data_in.append(1)
            row_in.append(n2)
            col_in.append(n1)
            nodes += [n1, n2]
        elif t == 1:
            in_edges.append(i)
            in_nodes2.append(n1)
            data.append(1)
            row.append(i)
            col.append(n1)
            data_in.append(1)
            row_in.append(n1)
            col_in.append(n2)
            data_vec_in.append(1)
            row_vec_in.append(i)
            col_vec_in.append(0)
        elif t == 2:
            out_edges.append(i)
            data.append(-1)
            row.append(i)
            col.append(n1)
            data_in.append(1)
            row_in.append(n1)
            col_in.append(n2)

    in_nodes2 = list(set(list(in_nodes2)))
    nodes = list(set(list(nodes)))
    for node in nodes:
        data_in.append(1)
        row_in.append(node)
        col_in.append(node)
    for edge in in_edges:
        for node in in_nodes:
            data.append(1)
            row.append(edge)
            col.append(node)
    for edge in out_edges:
        for node in out_nodes:
            data.append(1)
            row.append(edge)
            col.append(node)
    for node in in_nodes:
        data_diag.append(-1)
        row_diag.append(node)
        col_diag.append(node)
        for node2 in in_nodes2:
            data_in.append(1)
            row_in.append(node)
            col_in.append(node2)
    
    for node in out_nodes:
        data_diag.append(1)
        row_diag.append(node)
        col_diag.append(node)

    return spr.csr_matrix((data, (row, col)), shape=(ne, sid.nsq)), spr.csr_matrix((data_in, (row_in, col_in)), shape=(sid.nsq, sid.nsq)), spr.csr_matrix((data_diag, (row_diag, col_diag)), shape=(sid.nsq, sid.nsq)), spr.csc_matrix((data_vec_in, (row_vec_in, col_vec_in)), shape=(ne, 1))