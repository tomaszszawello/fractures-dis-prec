import numpy as np
import networkx as nx

def find_triangles(G):
    ''' Find triangles in the network (cells consisting of 3 nodes and 3 edges).
    Should keep their positions for plotting, maybe also edges.
    '''
    triangles = []
    triangles_pos = []
    pos = nx.get_node_attributes(G, 'pos')
    for edge in G.edges():
        n1, n2 = edge
        neigh1 = G.neighbors(n1)
        neigh2 = G.neighbors(n2)
        for node in set(neigh1).intersection(neigh2):
            triangles.append(tuple(sorted((n1, n2, node))))
    
    triangles = [*set(triangles)]
    for n1, n2, n3 in triangles:
        triangles_pos.append((np.array(pos[n1]) + np.array(pos[n2]) + np.array(pos[n3])) / 3)
    return triangles_pos