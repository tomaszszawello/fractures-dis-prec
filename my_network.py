import json
import networkx as nx
import numpy as np

from networkx.readwrite import json_graph

from config import simInputData

def load_my(sid:simInputData):
    fp = open(sid.load_name + '.json')
    G = json_graph.node_link_graph(json.load(fp))
    G2 = G.copy() # copy for saving network exactly same as initial, with only apertures changed
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
        
    # for n in G.nodes(): # hopefully it never happens
    #     if len(list(G.neighbors(n))) == 1:
    #         out_nodes.append(n)
    
    apertures = nx.get_edge_attributes(G, 'b').values()
    sid.b0 = sum(apertures) / len(apertures)
    lens = nx.get_edge_attributes(G, 'length').values()
    sid.l0 = sum(lens) / len(lens)
    sid.nsq = len(G.nodes())
    sid.n = int(np.sqrt(sid.nsq))
    return G, G2, in_nodes, out_nodes
