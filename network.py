""" Load network and manage all its properties.

This module contains classes and functions connected with loading the discrete
fracture network.

Notable classes
-------
Graph(nx.graph.Graph)
    container for network and its properties

Notable functions
-------
load(SimInputData) -> Graph
    load discrete fracture network from JSON file
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
    """
    in_nodes = []
    "list of inlet nodes"
    out_nodes = []
    "list of outlet nodes"

    def __init__(self):
        nx.graph.Graph.__init__(self)

    def update_network(self, edges: Edges) -> None:
        """ Updates apertures and flow in the graph.

        Parameters
        -------
        sid : SimInputData class object
            all config parameters of the simulation

        edges : Edges class object
            all edges in network and their parameters
        """
        nx.set_edge_attributes(self, dict(zip(edges.edge_list, \
            edges.apertures)), 'b')
        nx.set_edge_attributes(self, dict(zip(edges.edge_list, edges.flow)), \
            'q')

    def update_initial_network(self, sid: SimInputData, edges: Edges) -> None:
        """ Updates apertures and permeability in the initial graph.

        Parameters
        -------
        sid : SimInputData
            all config parameters of the simulation
        
        edges : Edges class object
            all edges in network and their parameters
        """
        nx.set_edge_attributes(self, dict(zip(edges.edge_list, \
            edges.apertures * sid.b0)), 'b')
        nx.set_edge_attributes(self, dict(zip(edges.edge_list, \
            (edges.apertures * sid.b0) ** 2 / 12)), 'perm')

    def dump_json_graph(self, sid: SimInputData, edges: Edges) -> None:
        """ Write graph out in json format.

        Parameters
        -------
        sid : SimInputData
            all config parameters of the simulation

        edges : Edges class object
            all edges in network and their parameters
        """
        # updata data in the graph
        self.update_initial_network(sid, edges)
        name = f'network_{sid.old_t:04f}'
        print("--> Dumping Graph into file: " + name + ".json")
        jsondata = json_graph.node_link_data(self)
        with open(sid.dirname + '/' + name + '.json', 'w') as fp:
            json.dump(jsondata, fp)
        print("--> Complete")


def load(sid:SimInputData) -> tuple[Graph, Graph]:
    """ Loads network from JSON file and translates it to Graph subclass.

    This function loads a discrete fracture network from file set in config and
    creates two instances of Graph subclass: graph and graph_real. graph_real
    is exactly like the loaded network - we only evolve its apertures and
    permeabilities. graph is slightly modified - nodes corresponding to source
    and target are removed and instead a group of inlet and outlet nodes
    (which in graph_real are connected to source and target by edges of
    infinite permeability) is created. graph is afterwards used to initialize
    all classes necessary for simulation.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation

    Returns
    -------
    graph : Graph class object
        network and all its properties

    graph_real : Graph class object
        network with structure exactly as loaded from file, but with evolved
        apertures and permeabilities
    """
    # load network from file and change it to Graph subclass
    fp = open(sid.load_name + '.json')
    graph = json_graph.node_link_graph(json.load(fp))
    graph.__class__ = Graph # TO DO: load into subclass more elegantly
    # copy for saving network exactly same as initial, but with evolving
    # apertures and permeabilities
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
    sid.n_nodes = len(graph.nodes())
    return graph, graph_real
