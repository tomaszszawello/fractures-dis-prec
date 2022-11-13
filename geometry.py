import delaunay as De

def set_geometry(n, geo='rect', *args, **kwargs):
    in_nodes, out_nodes= [], []
    if geo == 'rect':
        in_nodes = list(range(n))
        out_nodes = list(range(n * (n - 1), n * n))
    elif geo == 'own':
        in_nodes_pos = kwargs['in_nodes']
        out_nodes_pos = kwargs['out_nodes']
        for pos in in_nodes_pos:
            in_nodes.append(De.find_node(G, pos))
        for pos in out_nodes_pos:
            out_nodes.append(De.find_node(G, pos))
    return in_nodes, out_nodes
