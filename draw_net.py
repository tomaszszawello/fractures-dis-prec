import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import gridspec

from config import simInputData

def uniform_hist(sid:simInputData, G, in_nodes, out_nodes, boundary_edges, cb, cc, vols, triangles_pos, name):
    d_hist = []
    q_hist = []

    edges1 = []
    qs1 = []
    edges2 = []
    qs2 = []

    for n1, n2 in G.edges():
        q = G[n1][n2]['q']
        d = G[n1][n2]['d']

        q_hist.append(np.abs(q))
        d_hist.append(d)

        if (n1, n2) not in boundary_edges and (n2, n1) not in boundary_edges:
            if G[n1][n2]['c']:
                edges1.append((n1, n2))
                qs1.append(d)
            else:
                edges2.append((n1, n2))
                qs2.append(d)                

    pos = nx.get_node_attributes(G, 'pos')

    x_in, y_in = [], []
    for node in in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
        
    x_out, y_out = [], []
    for node in out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])

    # x_tr, y_tr = [], []
    # for node in triangles_pos:
    #     x_tr.append(node[0])
    #     y_tr.append(node[1])

    cols = 5
    plt.figure(figsize=(sid.figsize * 1.25, sid.figsize))
    spec = gridspec.GridSpec(ncols=cols, nrows=2, height_ratios=[5, 1])
    
    plt.subplot(spec.new_subplotspec((0, 0), colspan=cols))
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black', edgecolors='white')
    #plt.scatter(x_tr, y_tr, s=1*(vols < 9), facecolors='red', edgecolors='black')
    nx.draw_networkx_edges(G, pos, edgelist=edges1, edge_color = 'r', width=sid.ddrawconst * np.array(qs1))
    nx.draw_networkx_edges(G, pos, edgelist=edges2, edge_color = 'k', width=sid.ddrawconst * np.array(qs2))    
    #nx.draw_networkx_edges(G, pos, edgelist=edges, width=sid.ddrawconst * np.array(qs), edge_color=colors)
    #nx.draw_networkx_nodes(G, pos, node_size = 25 * oxdraw, node_color = oxdraw, cmap='Reds')
    plt.axis('equal')
    
    plt.subplot(spec[cols]).set_title('Diameter')
    plt.hist(d_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[cols + 1]).set_title('Flow')
    plt.hist(q_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[cols + 2]).set_title('cb')
    plt.hist(cb, bins=50)
    plt.yscale("log")

    plt.subplot(spec[cols + 3]).set_title('cc')
    plt.hist(cc, bins=50)
    plt.yscale("log")

    # plt.subplot(spec[cols + 4]).set_title('vola')
    # plt.hist(vols, bins=50)
    # plt.yscale("log")
    
    plt.savefig(sid.dirname + "/" + name)
    plt.close()

def uniform_hist(sid:simInputData, G, in_nodes, out_nodes, boundary_edges, cb, name, data):
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
            if data == 'd':
                qs.append(d)
            else:
                qs.append(sid.qdrawconst * q)                

    pos = nx.get_node_attributes(G, 'pos')

    x_in, y_in = [], []
    for node in in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
        
    x_out, y_out = [], []
    for node in out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])

    # x_tr, y_tr = [], []
    # for node in triangles_pos:
    #     x_tr.append(node[0])
    #     y_tr.append(node[1])

    cols = 5
    plt.figure(figsize=(sid.figsize * 1.25, sid.figsize))
    spec = gridspec.GridSpec(ncols=cols, nrows=2, height_ratios=[5, 1])
    
    plt.subplot(spec.new_subplotspec((0, 0), colspan=cols))
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black', edgecolors='white')
    #plt.scatter(x_tr, y_tr, s=1*(vols < 9), facecolors='red', edgecolors='black')
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color = 'k', width=sid.ddrawconst * np.array(qs))    
    #nx.draw_networkx_edges(G, pos, edgelist=edges, width=sid.ddrawconst * np.array(qs), edge_color=colors)
    nx.draw_networkx_nodes(G, pos, node_size = 25, node_color = 'k')
    plt.axis('equal')
    
    plt.subplot(spec[cols]).set_title('Diameter')
    plt.hist(d_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[cols + 1]).set_title('Flow')
    plt.hist(q_hist, bins=50)
    plt.yscale("log")

    plt.subplot(spec[cols + 2]).set_title('cb')
    plt.hist(cb, bins=50)
    plt.yscale("log")

    # plt.subplot(spec[cols + 3]).set_title('cc')
    # plt.hist(cc, bins=50)
    # plt.yscale("log")

    # plt.subplot(spec[cols + 4]).set_title('vola')
    # plt.hist(vols, bins=50)
    # plt.yscale("log")
    
    plt.savefig(sid.dirname + "/" + name)
    plt.close()