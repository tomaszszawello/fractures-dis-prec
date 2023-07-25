""" Initialize main simulation classes depending on config data.

This module creates instances of classes necessary for simulation. It loads 
a discrete fracture network and starts a new simulation with parameters set in
config file.

Notable functions
-------
build(SimInputData) -> tuple[Graph, Graph, In.Incidence, In.Edges, Data]
    create class objects and initialize their parameters
"""

import save as Sv

from config import SimInputData
from data import Data
import incidence as In
from network import Graph, load
from utils import make_dir


def build(sid:SimInputData) -> tuple[Graph, Graph, In.Incidence, In.Edges, \
    Data]:
    """ Initialize main classes used in simulation based on config file.

    Create class objects and initialize their parameters. Make a simulation
    directory and save there a config file.

    Parameters
    -------
    sid : SimInputData
        all config parameters of the simulation

    Returns
    -------
    graph : Graph class object
        network on which simulation is performed

    graph_real : Graph class object
        network with structure exactly as loaded from file, but with evolved
        apertures
    
    inc : Incidence class object
        matrices of incidence

    edges : Edges class object
        all edges in network and their parameters

    data : Data class object
        physical properties of the network measured during simulation
    """
    # load simulation network from file
    graph, graph_real = load(sid)
    # create simulation directory
    make_dir(sid)
    # translate network into incidence matrices and edge data vectors
    inc = In.Incidence()
    edges = In.create_matrices(sid, graph, inc)
    # initialize class for collecting data
    data = Data(sid, graph)
    # save file with initial parameters in simulation directory
    Sv.save_config(sid)
    return graph, graph_real, inc, edges, data