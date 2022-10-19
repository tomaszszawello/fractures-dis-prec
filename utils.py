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

def update_diameters(sid, edges, pnow, cb_now):
    breakthrough = False
    if sid.adaptive_dt:
        growth = np.zeros(len(edges))
        dt = [sid.dtmax]
        for i, e in enumerate(edges):
            n1, n2, d, l, t = e
            if pnow[n1] > pnow[n2]:  
                q = d ** 4 * (pnow[n1] - pnow[n2]) / l
                dd = q * cb_now[n1] / (sid.Da * l * d) * (1 - np.exp(-sid.Da / (1 + sid.G * d) * d * l / q))
            elif pnow[n2] > pnow[n1]:
                q = d ** 4 * (pnow[n2] - pnow[n1]) / l
                dd = q * cb_now[n2] / (sid.Da * l * d) * (1 - np.exp(-sid.Da / (1 + sid.G * d) * d * l / q))
            else:
                dd = 0
            growth[i] = dd
            if dd > 0:
                dt.append(sid.growth_rate * d / dd)
        growth *= np.min(dt)

        for i, e in enumerate(edges):
            n1, n2, d, l, t = e
            d += growth[i]
            if d < sid.dmin:
                d = sid.dmin
            if d > sid.dmax:
                d = sid.dmax
            edges[i] = (n1, n2, d, l, t)
            if t == 2 and d >= sid.dbreak:
                breakthrough = True
        return edges, np.min(dt), breakthrough
    
    else:
        dt = sid.dt
        for i, e in enumerate(edges):
            n1, n2, d, l, t = e
            if pnow[n1] > pnow[n2]:  
                q = d ** 4 * (pnow[n1] - pnow[n2]) / l
                dd = q * cb_now[n1] / (sid.Da * l * d) * (1 - np.exp(-sid.Da / (1 + sid.G * d) * d * l / q))
            elif pnow[n2] > pnow[n1]:
                q = d ** 4 * (pnow[n2] - pnow[n1]) / l
                dd = q * cb_now[n2] / (sid.Da * l * d) * (1 - np.exp(-sid.Da / (1 + sid.G * d) * d * l / q))
            else:
                dd = 0
            d += dt * dd
            if d < sid.dmin:
                d = sid.dmin
            if d > sid.dmax:
                d = sid.dmax
            edges[i] = (n1, n2, d, l, t)
            if t == 2 and d >= sid.dbreak:
                breakthrough = True
        return edges, dt, breakthrough

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

def collect_data(sid, edges, cb_now, pnow, iter):
    # data = []
    # for n1, n2, d, l, t in edges:
    #     data.append([cb_now[n1], cb_now[n2], d])
    # def save_data(name, data):
    #     success = 0
    #     while success != 1:
    #         try:
    #             f = open(sid.dirname+'/'+name+'.txt', 'a')
    #             np.savetxt(f, data)
    #             f.close()
    #             success = 1
    #         except PermissionError:
    #             pass
    success = 0
    while success != 1:
        try:
            f = open(sid.dirname+'/params.txt', 'a')
            f.write("Iter "+str(iter))
            f.write("\n")
            for n1, n2, d, l, t in edges:
                q = d ** 4 / l * np.abs(pnow[n1] - pnow[n2])
                Da_eff_l = sid.Da / (1 + sid.G * d) * d * l / q
                np.savetxt(f, [[cb_now[n1], cb_now[n2], d, q, l, Da_eff_l]], fmt = '%.6e')
                f.write("\n")
            f.write("\n")
            f.write("\n")
            f.close()
            success = 1
        except PermissionError:
            pass

    # save_data('params', data)
    # save_data('pin', pnow[in_nodes[0]])
    # save_data('oxnow', np.average(oxnow))
    # save_data('oxnow2', np.average(oxnow ** 2))
    # save_data('vnow', np.average(vnow))
    # save_data('vnow2', np.average(vnow ** 2))
    # save_data('S', S)
    # save_data('V', V)