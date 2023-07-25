import networkx as nx
import dill
import json

from networkx.readwrite import json_graph

from config import SimInputData
                                                                            
def save(name, sid:SimInputData, G, in_nodes, out_nodes, boundary_edges):
    pos = nx.get_node_attributes(G,'pos')
    All = [sid, [], in_nodes, out_nodes, boundary_edges, pos]

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
            G1.add_edge(n1, n2, d = d, q = 0, l = l)

        return G1
    
    G1 = reproduct()
    return  sid, G1, in_nodes, out_nodes, boundary_edges


def save_config(sid:SimInputData):
    f = open(sid.dirname+'/config.txt', 'w')
    for key, val in sid.__class__.__dict__.items():
        f.write(f'{key} = {val} \r')
    f.close()

def dump_json_graph(sid, G, name):
    """Write graph out in json format
 
    Parameters
    ---------- 
        self : object 
            DFN Class
        G :networkX graph
            NetworkX Graph based on the DFN
        name : string
             Name of output file (no .json)

    Returns
    -------

    Notes
    -----

"""
    print("--> Dumping Graph into file: " + name + ".json")
    jsondata = json_graph.node_link_data(G)
    with open(sid.dirname + '/' + name + '.json', 'w') as fp:
        json.dump(jsondata, fp)
    print("--> Complete")