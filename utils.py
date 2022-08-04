import networkx as nx
import numpy as np
import os
import scipy.sparse.linalg as sprlin
from delaunay import find_node


def PosGauss(mean, sigma):
    x = np.random.normal(mean, sigma)
    return(x if x>=0 else PosGauss(mean,sigma))
    
def mu_d(d):
    #return mu * (1 + 1 / (15 * d + 0.1))
    #return mu
    d = 100 * d + 0.2
    return ((6*np.exp(-0.085 * d)+2.2-2.44*np.exp(-0.06*d**0.645))*(d/(d-1.1)) ** 2 + 1) * (d/(d-1.1)) ** 2 / 1000

def solve_equation(matrix, result):
    return sprlin.spsolve(matrix, result)

def d_update(F, t):
    if (F > t.F0):
        if (F < t.F1):
            result = t.z0+(F-t.F0)*(t.z1-t.z0)/(t.F1-t.F0)
        else:
            result = t.z1
    else:
        result = t.z0
    return result

def update_diameters(sid, edges, d_pres, d_vegf, d_s):
    d_new = sid.c_pres * d_pres + sid.c_vegf * d_vegf + sid.c_s * d_s - sid.decrease
    for i,e in enumerate(edges):
        n1, n2, d, l, t = e
        if sid.linear:
            d += d_new[i] * d
        else:
            d += d_new[i]
        if d < sid.dmin:
            d = sid.dmin
        if d > sid.dmax:
            d = sid.dmax
        edges[i] = (n1, n2, d, l, t)
    return edges

def make_dir(sid):
        # if not os.path.isdir(sid.dirname):
        #     os.makedirs(sid.dirname)

        i = 0
        dirname2 = sid.dirname
        while (sid.dirname == dirname2):
            if not os.path.isdir(sid.dirname + "/" + str(i)):
                sid.dirname = sid.dirname + "/" + str(i)
            else:
                i += 1
        os.makedirs(sid.dirname)

class fParams():
    def __init__(self, params):
        self.F0 = params['F0']
        self.F1 = params['F1']
        self.z0 = params['z0']
        self.z1 = params['z1']


def find_ladder(sid, G):
    # in1 = find_node(G, [40, 40])
    # in2 = find_node(G, [25, 65])
    # out1 = find_node(G, [60, 60])
    # out2 = find_node(G, [35, 75])
    #for i in range(len(p1) - 1):
    #    G[p1[i]][p1[i+1]]['d'] = sid.dth
    #for i in range(len(p2) - 1):
    #    G[p2[i]][p2[i+1]]['d'] = sid.dth
    in1 = find_node(G, [40, 40])
    in2 = find_node(G, [30, 50])
    in3 = find_node(G, [20, 60])
    out1 = find_node(G, [60, 60])
    out2 = find_node(G, [50, 70])
    out3 = find_node(G, [40, 80])
    p1 = nx.shortest_path(G, in1, in2)
    p2 = nx.shortest_path(G, out1, out2)
    p3 = nx.shortest_path(G, in2, in3)
    p4 = nx.shortest_path(G, out2, out3)
    p5 = nx.shortest_path(G, in2, out2)
    p6 = nx.shortest_path(G, in3, out3)
    for i in range(len(p1) - 1):
       G[p1[i]][p1[i+1]]['d'] = sid.dth
    for i in range(len(p2) - 1):
       G[p2[i]][p2[i+1]]['d'] = sid.dth
    for i in range(len(p3) - 1):
       G[p3[i]][p3[i+1]]['d'] = sid.dth
    for i in range(len(p4) - 1):
       G[p4[i]][p4[i+1]]['d'] = sid.dth
    for i in range(len(p5) - 1):
       G[p5[i]][p5[i+1]]['d'] = sid.dth
    for i in range(len(p6) - 1):
       G[p6[i]][p6[i+1]]['d'] = sid.dth
    in_nodes = [in1]
    out_nodes = [out1]
    in_nodes_ox = [in1]
    out_nodes_ox = [out1]
    oxresult = np.zeros(sid.nsq)
    for node in p1 + p2 + p3 + p4 + p5 + p6:
        oxresult[node] = 1
    return G, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, oxresult

class simAnalysisData:
    pressure = []
    oxygen = []
    vegf = []
    signal = []

def collect_data(sid, edges, in_nodes, out_nodes, pnow, vnow, oxnow, oxresult):
    V = 0
    V_q = 0
    S = 0
    S_q = 0
    q2 = 0
    q4 = 0
    N = 0
    L = 0
    for n1, n2, d, l, t in edges:
        q = sid.c1 / sid.mu * d ** 4 * np.abs(pnow[n1] - pnow[n2]) / l
        V0 = np.pi * (d / 2) ** 2 * l
        S0 = np.pi * d * l
        V += V0
        S += S0
        if q > sid.q_prun:
            V_q += V0
            S_q += S0
            q2 += q ** 2
            q4 += q ** 4
            N += 1
        if (oxresult[n1] == 1 and oxresult[n2] == 1) or (oxresult[n1] == 2 and oxresult[n2] == 2):
            L += l
    A = (N - q2 ** 2 / q4) / (N - 1)
    ox_out = 0
    for node in out_nodes:
        ox_out += oxnow[node]
    ox_out = ox_out / len(out_nodes)
    oxnow_nodes1 = np.sum(oxnow>0.1)
    oxnow_nodes2 = np.sum(oxnow>0.2)
    oxnow_nodes3 = np.sum(oxnow>0.3)
    oxnow_nodes4 = np.sum(oxnow>0.4)
    oxnow_nodes5 = np.sum(oxnow>0.5)
    oxnow_nodes6 = np.sum(oxnow>0.6)
    oxnow_nodes7 = np.sum(oxnow>0.7)
    oxnow_nodes8 = np.sum(oxnow>0.8)
    oxnow_nodes9 = np.sum(oxnow>0.9)
    oxresult_nodes = np.sum(oxresult>0)
    data = [pnow[in_nodes[0]], np.average(oxnow), np.average(oxnow ** 2), np.average(vnow), np.average(vnow ** 2), V, S, V_q, ox_out, A, oxnow_nodes1, oxnow_nodes2, oxnow_nodes3, oxnow_nodes4, oxnow_nodes5, oxnow_nodes6, oxnow_nodes7, oxnow_nodes8, oxnow_nodes9, oxresult_nodes, L]
    def save_data(name, data):
        success = 0
        while success != 1:
            try:
                f = open(sid.dirname+'/'+name+'.txt', 'a')
                np.savetxt(f, [data])
                f.close()
                success = 1
            except PermissionError:
                pass
    save_data('params', data)
    # save_data('pin', pnow[in_nodes[0]])
    # save_data('oxnow', np.average(oxnow))
    # save_data('oxnow2', np.average(oxnow ** 2))
    # save_data('vnow', np.average(vnow))
    # save_data('vnow2', np.average(vnow ** 2))
    # save_data('S', S)
    # save_data('V', V)