import networkx as nx
import dill
from config import simInputData
from utils import fParams
                                                                            
def save(name, sid:simInputData, G, edges, in_nodes, out_nodes, boundary_edges):
    pos = nx.get_node_attributes(G,'pos')
    All = [sid, edges, in_nodes, out_nodes, boundary_edges, pos]

    with open(sid.dirname+name, 'wb') as file:
        dill.dump(All, file)


def load(name):
    with open(name, 'rb') as file:
        All= dill.load(file)
    
    sid = All[0]
    edges = All[1]
    in_nodes = All[2]
    out_nodes = All[3]
    boundary_edges = All[4]
    pos = All[5]

    def reproduct():
        G1 = nx.Graph()
        new_pos = {}
        for key, value in pos.items():
            new_pos[int(key)] =  value
        
        for node in new_pos:
            G1.add_node(node, pos = new_pos[node])
        for n1, n2, d, l, t in edges:
            G1.add_edge(n1, n2, d = d, q = 0, length = l)

        return G1
    
    G1 = reproduct()
    return  sid, G1, edges, in_nodes, out_nodes, boundary_edges


def save_config(sid:simInputData):
    f = open(sid.dirname+'/config.txt', 'w')
    for key, val in sid.__class__.__dict__.items():
        if isinstance(val, fParams):
            f.write(f'{key} = {val.__dict__} \r')
        else:
            f.write(f'{key} = {val} \r')
    f.close()     