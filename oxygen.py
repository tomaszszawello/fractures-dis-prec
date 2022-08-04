import numpy as np
import scipy.sparse as spr

from utils import d_update

from config import simInputData


def create_vector(sid:simInputData, in_nodes_ox, out_nodes_ox):
    oxresult = np.zeros(sid.nsq)
    for node in in_nodes_ox:
        oxresult[node] = 1
    for node in out_nodes_ox:
        oxresult[node] = 2
    return oxresult

def create_blood_vector(sid:simInputData, in_nodes_ox):
    bloodoxresult = np.zeros(sid.nsq)
    for node in in_nodes_ox:
        bloodoxresult[node] = sid.cox_in
    return bloodoxresult 

def update_matrix(sid:simInputData, oxresult, bloodoxresult, pnow, edges):
    data, row, col = [], [], []

    diag = np.ones(sid.nsq)
    for n1, n2, d, l, t in edges:
        if bloodoxresult[n1] == 1:
            diag[n1] = 1
        elif oxresult[n1] != 0:
            diag[n1] = np.pi * d * sid.k_ox_b
        else:
            diag[n1] = sid.k_ox_t
        if bloodoxresult[n2] == 1:
            diag[n2] = 1
        elif oxresult[n2] != 0:
            diag[n2] = np.pi * d * sid.k_ox_b
        else:
            diag[n2] = sid.k_ox_t         
    for n1, n2, d, l, t in edges:
        if bloodoxresult[n1] == 1:
            if oxresult[n2] == 1:
                if pnow[n2] < pnow[n1]:
                    res = -(pnow[n1] - pnow[n2]) * d ** 4 * np.pi / (128 * sid.mu * l) #
                    data.append(res)
                    row.append(n2)
                    col.append(n1)
                    diag[n2] -= res
                res = -np.pi * d ** 2 / 4 * sid.D_ox_b / l
            else:
                res = -sid.D_ox_t / l
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n2] -= res
        elif bloodoxresult[n2] == 1:
            if oxresult[n1] == 1:
                if pnow[n1] < pnow[n2]:
                    res = -(pnow[n2] - pnow[n1]) * d ** 4 * np.pi / (128 * sid.mu * l)
                    data.append(res)
                    row.append(n1)
                    col.append(n2)
                    diag[n1] -= res
                res = -np.pi * d ** 2 / 4 * sid.D_ox_b / l
            else:
                res = -sid.D_ox_t / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
        elif oxresult[n1] == 1 and oxresult[n2] == 1:
            if pnow[n1] < pnow[n2]:
                res = -(pnow[n2] - pnow[n1]) * d ** 4 * np.pi / (128 * sid.mu * l)
                data.append(res)
                row.append(n1)
                col.append(n2)
                diag[n1] -= res
            elif pnow[n2] < pnow[n1]:
                res = -(pnow[n1] - pnow[n2]) * d ** 4 * np.pi / (128 * sid.mu * l)
                data.append(res)
                row.append(n2)
                col.append(n1)
                diag[n2] -= res
            res = -np.pi * d ** 2 / 4 * sid.D_ox_b / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n2] -= res
        elif oxresult[n1] == 2 and oxresult[n2] == 2:
            if pnow[n1] < pnow[n2]:
                res = -(pnow[n2] - pnow[n1]) * d ** 4 * np.pi / (128 * sid.mu * l)
                data.append(res)
                row.append(n1)
                col.append(n2)
                diag[n1] -= res
            elif pnow[n2] < pnow[n1]:
                res = -(pnow[n1] - pnow[n2]) * d ** 4 * np.pi / (128 * sid.mu * l)
                data.append(res)
                row.append(n2)
                col.append(n1)
                diag[n2] -= res
            res = -np.pi * d ** 2 / 4 * sid.D_ox_b / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n2] -= res
        else:
            res = -sid.D_ox_t / l
            data.append(res)
            row.append(n1)
            col.append(n2)
            diag[n1] -= res
            data.append(res)
            row.append(n2)
            col.append(n1)
            diag[n2] -= res

    for node, datum in enumerate(diag):
        if datum != 0:
            row.append(node)
            col.append(node)
            data.append(datum)

    return spr.csr_matrix((data, (row, col)), shape=(sid.nsq, sid.nsq))


def update_graph(sid:simInputData, oxnow, oxresult, edges):
    oxresult2 = oxresult.copy()
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        if (oxresult2[n1] == 1 or oxresult2[n2] == 1):
            F = sid.F_mult_ox * np.abs(oxnow[n1] - oxnow[n2])/l
            d += d_update(F, sid.F_ox)
            if d > sid.dth:
                oxresult[n1] = 1
                oxresult[n2] = 1
        edges[i] = (n1, n2, d, l)


    return edges, oxresult

def update_oxresult(sid:simInputData, edges, oxresult):
    for e in edges:
        n1, n2, d, l, t = e
        if (oxresult[n1] == 1 and oxresult[n2] == 1):
            if d < sid.dth:
                oxresult[n1] = 0
                oxresult[n2] = 0
        elif (oxresult[n1] == 1 or oxresult[n2] == 1):
            if d > sid.dth:
                oxresult[n1] = 1
                oxresult[n2] = 1
        
    return oxresult