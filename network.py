""" Build network and manage all its properties.

This module contains classes and functions connected with building Delaunay
network, setting boundary condition on it and evolving.

Notable classes
-------
Graph(nx.graph.Graph)
    container for network and its properties

Notable functions
-------
load(SimInputData) -> Graph
    build Delaunay network with parameters from config
"""

from __future__ import annotations
from networkx.readwrite import json_graph
from typing import TYPE_CHECKING
import json
import networkx as nx
import numpy as np

from config import SimInputData
if TYPE_CHECKING:
    from incidence import Edges


class Graph(nx.graph.Graph):
    """ Contains network and all its properties.

    This class is derived from networkx Graph and contains all information
    abount the network and its properties.

    Attributes
    -------
    in_nodes : list
        list of inlet nodes
    out_nodes : list
        list of outlet nodes
    """
    in_nodes = []
    out_nodes = []

    def __init__(self):
        nx.graph.Graph.__init__(self)

    def update_network(self, edges: Edges) -> None:
        """ Update diameters and flow in the graph.

        Parameters
        -------
        edges : Edges class object
            all edges in network and their parameters
            edge_list - array of tuples (n1, n2) with n1, n2 being nodes
            connected by edge with a given index
            diams - diameters of edges
            flow - flow in edges
        """
        nx.set_edge_attributes(self, dict(zip(edges.edge_list, \
            edges.apertures)), 'b')
        nx.set_edge_attributes(self, dict(zip(edges.edge_list, edges.flow)), \
            'q')

    def update_initial_network(self, sid: SimInputData, edges: Edges) -> None:
        """ Update diameters and flow in the graph.

        Parameters
        -------
        edges : Edges class object
            all edges in network and their parameters
            edge_list - array of tuples (n1, n2) with n1, n2 being nodes
            connected by edge with a given index
            diams - diameters of edges
            flow - flow in edges
        """
        nx.set_edge_attributes(self, dict(zip(edges.edge_list, \
            edges.apertures * sid.b0)), 'b')
        nx.set_edge_attributes(self, dict(zip(edges.edge_list, \
            (edges.apertures * sid.b0) ** 2 / 12)), 'perm')


def load(sid:SimInputData) -> tuple[Graph, Graph]:
    # load network from file and change it to Graph subclass
    fp = open(sid.load_name + '.json')
    graph = json_graph.node_link_graph(json.load(fp))
    graph.__class__ = Graph # TO DO: load into subclass more elegantly
    # copy for saving network exactly same as initial, buut with evolving
    # apertures
    graph_real = graph.copy()
    # get rid of additional edges with inf permeability from inlet/outlet node
    # (we keep them in graph_real)
    for edge in graph.edges():
        n1, n2 = edge
        if isinstance(n1, str) or isinstance(n2, str):
            if n1 == 's':
                graph.in_nodes.append(n2)
            elif n1 == 't':
                graph.out_nodes.append(n2)
            if n2 == 's':
                graph.in_nodes.append(n1)
            elif n2 == 't':
                graph.out_nodes.append(n1)
            graph.remove_edge(n1, n2)
    # get rid of additional inlet/outlet node
    remove_nodes = []
    for node in graph.nodes():
        if isinstance(node, str):
            remove_nodes.append(node)
    for node in remove_nodes:
        graph.remove_node(node)
    # update parameters in sid based on loaded graph
    apertures = nx.get_edge_attributes(graph, 'b').values()
    sid.b0 = sum(apertures) / len(apertures)
    lens = nx.get_edge_attributes(graph, 'length').values()
    nx.set_edge_attributes(graph, 0, 'q')
    sid.l0 = sum(lens) / len(lens)
    sid.nsq = len(graph.nodes())
    sid.n = int(np.sqrt(sid.nsq))
    return graph, graph_real
