""" Perform particle tracking.

This module contains function for particle tracking - standard, weighted by
the concentration of reactant or with removing particles due to reaction.

Notable functions
-------
track(SimInputData, Incidence, Graph, Edges, Data, numpy ndarray) -> None
    perform tracking and collect relative data
"""

import numpy as np
import scipy.sparse as spr

from config import SimInputData
from data import Data
from network import Graph
from incidence import Edges, Incidence

def track(sid: SimInputData, inc: Incidence, graph: Graph, \
        edges: Edges, data: Data, pressure: np.ndarray) -> None:
    """ Perform particle tracking and collect flow/velocity data.

    This function performs particle tracking and saves them in Data class to
    later create breakthrough curves. Depending on config parameters, it
    performs standard tracking, concentration weighted tracking and tracking
    with removing particles due to reactions (where in each edge, we remove
    a tracked particle with probability dependent on reaction term in a given
    edge - we calculate how much the concentration changes in this edge
    c_out / c_in = exp(- 2 Da / (1 + G b) L / q)) and we remove the particle
    with probability p = 1 - c_out / c_in). We also collect the flow and
    velocity in the whole network.
    
    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation

    inc : Incidence class object
        matrices of incidence

    graph : Graph class object
        network and all its properties

    edges : Edges class object
        all edges in network and their parameters
    
    data : Data class object
        physical properties of the network measured during simulation
        
    pressure : numpy ndarray
        vector of pressure in nodes
    """
    data.tracking_names.append(sid.old_t)
    breakthrough_times = []
    concentrations = []
    # find upstream neighbours
    neigh_inc = (spr.diags(edges.flow) @ inc.incidence > 0).T
    tot_flow = np.abs(edges.flow * edges.fracture_lens)
    tot_velocity = np.abs(edges.flow / edges.apertures)
    # collect data for flow and velocity
    data.velocities.append(tot_velocity)
    data.vol_flow.append(tot_flow)
    # calculate travel time through each edge
    tot_time = np.abs(edges.lens / tot_velocity)
    # we introduce a particle to an inlet edge with probability proportional to
    # the flow in that edge
    inlet_flow = edges.fracture_lens * edges.apertures ** 3 \
        / edges.lens * (inc.inlet @ pressure)
    inlet_flow /= np.sum(inlet_flow)
    # standard and concentration weighted tracking
    if sid.include_tracking:
        # reaction term for calculation of concentration drop during tracking
        exp = np.exp(-np.abs(2 * sid.Da / (1 + sid.G * edges.apertures) \
            * edges.lens / edges.flow))
        # loop for tracking particles
        for _ in range(sid.n_part):
            time = 0
            conc = 1
            # choose inlet edge to introduce particle
            in_edge = np.random.choice(len(inlet_flow), p = inlet_flow)
            n1, n2 = edges.edge_list[in_edge]
            # put particle on the end of inlet edge with lower pressure
            if pressure[n1] > pressure[n2]:
                node = n2
            else:
                node = n1
            # travel until particle reaches an outlet node
            while node not in graph.out_nodes:
                prob = []
                # get edge indices to neighbors of nodes
                neigh_edges = neigh_inc[node].nonzero()[1]
                for edge in neigh_edges:
                    prob.append(tot_flow[edge])
                prob = np.array(prob) / np.sum(prob)
                # choose neighbor with probability dependent on flow
                edge = neigh_edges[np.random.choice(len(prob), p = prob)]
                # increase time and decrease concentration
                time += tot_time[edge]
                conc *= exp[edge]
                # if concentration is too low, reduce it to 0 (for plotting)
                if conc < 1e-40:
                    conc = 0
                # change node to the chosen one
                n1, n2 = edges.edge_list[edge]
                if n1 == node:
                    node = n2
                else:
                    node = n1
            breakthrough_times.append(time)
            concentrations.append(conc)
        # collect tracking times
        data.tracking.append(breakthrough_times)
        data.concentration_tracking.append(concentrations)
    # reactive tracking
    if sid.include_reactive_tracking:
        # scale exponent so as not to kill particles too fast
        exp = np.exp(-np.abs(2 * sid.Da / (1 + sid.G * edges.apertures) \
            * edges.lens / edges.flow) / 10)
        reactive_breakthrough_times = []
        # loop for tracking particles
        for _ in range(sid.n_part):
            time = 0
            conc = 1
            # choose inlet edge to introduce particle
            in_edge = np.random.choice(len(inlet_flow), p = inlet_flow)
            n1, n2 = edges.edge_list[in_edge]
            # put particle on the end of inlet edge with lower pressure
            if pressure[n1] > pressure[n2]:
                node = n2
            else:
                node = n1
            # travel until particle is killed or reaches an outlet node
            flag = True
            while flag:
                prob = []
                # get edge indices to neighbors of nodes
                neigh_edges = neigh_inc[node].nonzero()[1]
                for edge in neigh_edges:
                    prob.append(tot_flow[edge])
                prob = np.array(prob) / np.sum(prob)
                # choose neighbor with probability dependent on flow
                edge = neigh_edges[np.random.choice(len(prob), p = prob)]
                time += tot_time[edge]
                conc *= exp[edge]
                # kill particle with probability depending on amount of
                # reaction in a given edge
                if np.random.rand() < conc * (1 - exp[edge]):
                    # if particle is killed, break loop and skip it in data
                    break
                n1, n2 = edges.edge_list[edge]
                # change node to the chosen one
                if n1 == node:
                    node = n2
                else:
                    node = n1
                # if particle reached outlet, end loop
                if node in graph.out_nodes:
                    flag = False
            # if particle reached outlet, include it in data
            if not flag:
                reactive_breakthrough_times.append(time)
        # collect reactive tracking times
        data.reactive_tracking.append(reactive_breakthrough_times)
