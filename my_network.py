import json
import networkx as nx
import numpy as np

from networkx.readwrite import json_graph

from config import simInputData

# def build_my(sid:simInputData):
#     n = 6
#     sid.nsq = n
#     G = nx.Graph()
#     G.add_nodes_from(list(range(n)))

#     points = [(3, 3), (2.5, 1.5), (1.75, 2.25), (4, 4), (1.25, 0), (0.75, 1)]
#     edges = [(0, 3), (0, 1), (0, 2), (1, 2), (1, 4), (2, 5)]
#     diams = [3.4641016151377547e-06, 3.4641016151377547e-06, 3.4641016151377547e-06, 3.4641016151377547e-06, 4.898979485566356e-06, 3.4641016151377547e-06]
#     lens = [0.5, 0.2, 0.2, 0.4, 0.5, 0.5]
#     fl = np.ones(n)

#     sid.d0 = np.average(diams)
#     sid.l0 = np.average(lens)
#     diams = diams / sid.d0
#     lens = lens / sid.l0

#     for node in G.nodes:
#         G.nodes[node]['pos'] = points[node]
#         G.nodes[node]['fl'] = fl[node]

#     G.add_edges_from(edges)

#     for i, e in enumerate(G.edges()):
#         n1, n2 = e
#         G[n1][n2]['d'] = diams[i]
#         G[n1][n2]['l'] = lens[i]

#     in_nodes = [3]
#     out_nodes = [4, 5]
#     boundary_edges = []
#     return G, in_nodes, out_nodes, boundary_edges

def load_my(sid:simInputData):
    fp = open(sid.load_name + '.json')
    G = json_graph.node_link_graph(json.load(fp))
    in_nodes, out_nodes = [], []
    for e in G.edges():
        n1, n2 = e
        if isinstance(n1, str) or isinstance(n2, str):
            if n1 == 's':
                in_nodes.append(n2)
            elif n1 == 't':
                out_nodes.append(n2)
            if n2 == 's':
                in_nodes.append(n1)
            elif n2 == 't':
                out_nodes.append(n1)
            G.remove_edge(n1, n2)
    remove_nodes = []
    for n in G.nodes():
        if isinstance(n, str):
            remove_nodes.append(n)
    for n in remove_nodes:
        G.remove_node(n)
    # for n in G.nodes():
    #     if len(list(G.neighbors(n))) == 1:
    #         out_nodes.append(n)
    apertures = nx.get_edge_attributes(G, 'b').values()
    sid.b0 = sum(apertures) / len(apertures)
    lens = nx.get_edge_attributes(G, 'length').values()
    sid.l0 = sum(lens) / len(lens)
    sid.nsq = len(G.nodes())
    sid.n = int(np.sqrt(sid.nsq))
    return G, in_nodes, out_nodes
