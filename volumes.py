import numpy as np
import networkx as nx
import scipy.sparse as spr

from collections import defaultdict

from config import simInputData

def find_triangles(G):
    ''' Find triangles in the network (cells consisting of 3 nodes and 3 edges).
    Should keep their positions for plotting, maybe also edges.
    '''
    triangles_dict = defaultdict(list)
    triangles = []
    triangles_pos = []
    triangles_inc_row = []
    triangles_inc_col = []
    triangles_inc_data = []
    pos = nx.get_node_attributes(G, 'pos')
    ne = len(G.edges())
    for i, e in enumerate(G.edges()):
        n1, n2 = e
        neigh1 = G.neighbors(n1)
        neigh2 = G.neighbors(n2)
        for node in set(neigh1).intersection(neigh2):
            triangles_dict[(tuple(sorted((n1, n2, node))))].append(i)
    nt = len(triangles_dict.keys())
    #triangles = [*set(triangles)]
    for i, nodes in enumerate(triangles_dict.keys()):
        n1, n2, n3 = nodes
        triangles_pos.append((np.array(pos[n1]) + np.array(pos[n2]) + np.array(pos[n3])) / 3)
        for edge in triangles_dict[nodes]:
            triangles_inc_row.append(edge)
            triangles_inc_col.append(i)
            triangles_inc_data.append(1)


    return spr.csr_matrix((triangles_inc_data, (triangles_inc_row, triangles_inc_col)), shape=(ne, nt)), triangles_pos

def create_vector(sid:simInputData, triangle_pos):
    ''' Create vector with initial volumes of A.
    '''
    return sid.vol_a_in * np.ones(len(triangle_pos))