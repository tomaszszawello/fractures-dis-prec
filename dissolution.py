import numpy as np
import scipy.sparse as spr

from config import simInputData


def create_vector(sid:simInputData, in_nodes):
    cb_result = np.zeros(sid.nsq)
    for node in in_nodes:
        cb_result[node] = sid.cb_in
    return cb_result

def update_matrix(sid:simInputData, pnow, edges):
    data, row, col = [], [], []
    diag = np.zeros(sid.nsq)
    for n1, n2, d, l, t in edges:
        if t == 0:
            if pnow[n1] > pnow[n2]:
                q = (pnow[n1] - pnow[n2]) * d ** 4 / l
                qc = q * np.exp(-sid.Da / (1 + sid.G * d) * d * l / q)
                data.append(qc)
                row.append(n2)
                col.append(n1)
                diag[n2] -= q
            elif pnow[n2] > pnow[n1]:
                q = (pnow[n2] - pnow[n1]) * d ** 4 / l
                qc = q * np.exp(-sid.Da / (1 + sid.G * d) * d * l / q)
                data.append(qc)
                row.append(n1)
                col.append(n2)
                diag[n1] -= q
        elif t == 1:
            q = np.abs(pnow[n2] - pnow[n1]) * d ** 4 / l
            qc = q * np.exp(-sid.Da / (1 + sid.G * d) * d * l / q)
            data.append(qc)
            row.append(n1)
            col.append(n2)
            diag[n1] -= q
            diag[n2] = 1
        elif t == 2:
            q = np.abs(pnow[n1] - pnow[n2]) * d ** 4 / l
            qc = q * np.exp(-sid.Da / (1 + sid.G * d) * d * l / q)
            data.append(qc)
            row.append(n2)
            col.append(n1)
            diag[n2] -= q
            
    for node, datum in enumerate(diag):
        if datum != 0:
            row.append(node)
            col.append(node)
            data.append(datum)
        else:
            row.append(node)
            col.append(node)
            data.append(1)
    return spr.csr_matrix((data, (row, col)), shape=(sid.nsq, sid.nsq))