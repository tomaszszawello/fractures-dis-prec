import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import gridspec

from config import simInputData

def uniform_hist(sid:simInputData, G, in_nodes, out_nodes, boundary_edges, cb, cc, name):
    d_hist = []
    q_hist = []

    edges = []
    qs = []

    for n1, n2 in G.edges():
        q = G[n1][n2]['q']
        d = G[n1][n2]['d']

        q_hist.append(np.abs(q))
        d_hist.append(d)

        if (n1, n2) not in boundary_edges and (n2, n1) not in boundary_edges:
            edges.append((n1, n2))
            qs.append(d)


    pos = nx.get_node_attributes(G, 'pos')

    x_in, y_in = [], []
    for node in in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
        
    x_out, y_out = [], []
    for node in out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])

    plt.figure(figsize=(sid.figsize * 1.25, sid.figsize))
    spec = gridspec.GridSpec(ncols=4, nrows=2, height_ratios=[5, 1])
    
    plt.subplot(spec.new_subplotspec((0, 0), colspan=4))
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black', edgecolors='white')
    #nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color = colors, width=sid.ddrawconst * np.array(qs), edge_cmap = plt.cm.plasma)
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=sid.ddrawconst * np.array(qs))
    #nx.draw_networkx_nodes(G, pos, node_size = 25 * oxdraw, node_color = oxdraw, cmap='Reds')   
    plt.axis('equal')
    
    plt.subplot(spec[4]).set_title('Diameter')
    plt.hist(d_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[5]).set_title('Flow')
    plt.hist(q_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[6]).set_title('cb')
    plt.hist(cb, bins=50)
    plt.yscale("log")

    plt.subplot(spec[7]).set_title('cc')
    plt.hist(cc, bins=50)
    plt.yscale("log")
    
    plt.savefig(sid.dirname + "/" + name)
    plt.close()
