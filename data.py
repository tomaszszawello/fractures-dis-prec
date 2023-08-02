""" Collect physical data from the simulation and save/plot them.

This module contains Data class, storing information about physical data in
the simulation. It stores the data during simulation and afterwards saves them
in a text file and plots them. For now the data are: permeability and
channelization (take a slice of the system in a given x-coordinate and check
how many of the edges contain half of the total flow through the slice) for
3 slices of the system (at 1/4, 1/2 and 3/4).

Notable classes
-------
Data
    container for physical data collected during simulation

TO DO: name data on plots, fix no values on y-axis for permeability in
channelization plot, maybe add plots of effluent concentration and dissolved
volume, include tracking plots in data?
"""

from matplotlib import gridspec
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.sparse as spr

from config import SimInputData
from incidence import Edges, Incidence
from network import Graph


class Data():
    """ Contains data collected during the simulation.

    This class is a container for all data collected during simulation, such
    as permeability, channelization etc. Part of the data
    is saved in params.txt file, the rest is plotted and saved as figures.
    """
    t: list = []
    "elapsed time of the simulation"
    perm: list = []
    "permeability between inlet and outlet"
    order: list = []
    "order parameter ((N - sum(q^2)^2/sum(q^4)) / (N - 1))"
    channels_1: list = []
    "channelization in slice 1 (defaultly 1/4 of system)"
    channels_2: list = []
    "channelization in slice 2 (defaultly 1/2 of system)"
    channels_3: list = []
    "channelization in slice 3 (defaultly 3/4 of system)"

    def __init__(self, sid: SimInputData, graph: Graph):
        self.dirname = sid.dirname
        # set positions of 3 default slices measured vs time
        pos_x = np.array(list(nx.get_node_attributes(graph, 'x').values()), \
            dtype = float)
        self.slice_x1 = float((np.min(pos_x) + np.average(pos_x)) / 2)
        self.slice_x2 = float(np.average(pos_x))
        self.slice_x3 = float((np.max(pos_x) + np.average(pos_x)) / 2)

    def collect(self, sid: SimInputData, graph: Graph, inc: Incidence, \
        edges: Edges, pressure: np.ndarray) -> None:
        """ Collect data from different vectors.

        This function extracts information such as permeability, channelization
        etc. and saves them in the data class.

        Parameters
        -------
        sid : SimInputData class object
            all config parameters of the simulation

        inc : Incidence class object
            matrices of incidence

        edges : Edges class object
            all edges in network and their parameters

        graph : Graph class object
            network and all its properties

        pressure : numpy ndarray
            vector of current pressure
        """
        # simulation time
        self.t.append(sid.old_t)
        # permeability (dimensionless)
        self.perm.append(1 / np.max(pressure))
        # flow focusing parameter
        self.order.append((np.sum(edges.flow != 0) \
            - np.sum(np.abs(edges.fracture_lens * edges.flow) ** 2) ** 2 \
            / np.sum(np.abs(edges.fracture_lens * edges.flow) ** 4)) \
            / (np.sum(edges.flow != 0) - 1))
        # channelization
        self.channels_1.append(self.check_channelization(graph, inc, edges, \
            self.slice_x1)[1])
        self.channels_2.append(self.check_channelization(graph, inc, edges, \
            self.slice_x2)[1])
        self.channels_3.append(self.check_channelization(graph, inc, edges, \
            self.slice_x3)[1])
        
    def plot(self) -> None:
        """ Save all data and plot them.        
        """
        self.save()
        self.plot_params()
        self.plot_channelization()

    def save(self) -> None:
        """ Saves data to text file.

        This function saves the collected data to text file params.txt in
        columns. If the simulation is continued from saved parameters, new data
        is appended to that previously collected.
        """
        is_saved = False
        while not is_saved: # prevents problems with opening text file
            try:
                file = open(self.dirname + '/params.txt', 'a', \
                    encoding = "utf-8")
                np.savetxt(file, np.array([self.t, self.perm, self.order,
                    self.channels_1, self.channels_2, self.channels_3], \
                    dtype = float).T)
                file.close()
                is_saved = True
            except PermissionError:
                pass

    def check(self, edges: Edges) -> None:
        """ Check if flow in the system is valid.

        This function calculates and prints inflow and outflow to check if they
        are equal.

        Parameters
        -------
        edges : Edges class object
            all edges in network and their parameters
        """
        Q_in = np.sum(edges.inlet * np.abs(edges.fracture_lens * edges.flow))
        Q_out = np.sum(edges.outlet * np.abs(edges.fracture_lens * edges.flow))
        print('Q_in =', Q_in, 'Q_out =', Q_out)

    def plot_params(self) -> None:
        """ Plots data from text file.

        This function loads the data from text file params.txt and plots them
        to file params.png.
        """
        f = open(self.dirname + '/params.txt', 'r', encoding = "utf-8")
        data = np.loadtxt(f)
        n_data = data.shape[1]
        # first column is time
        t = data[:, 0]
        plt.figure(figsize = (15, 5))
        plt.suptitle('Parameters')
        spec = gridspec.GridSpec(ncols = n_data - 1, nrows = 1)
        for i_data in range(n_data - 1):
            # TO DO: name data columns
            plt.subplot(spec[i_data]).set_title(f'Data {i_data}')
            plt.plot(t, data[:, i_data + 1])
            plt.xlabel('simulation time')
        plt.savefig(self.dirname + '/params.png')
        plt.close()
        
    def plot_channelization(self) -> None:
        """ Plot channelization data from text file.

        This function loads the data from text file params.txt and plots those
        corresponding to channelization vs permeability.
        """
        f = open(self.dirname+'/params.txt', 'r')
        data = np.loadtxt(f)
        n_data = data.shape[1]
        t = data[:, 0]
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('simulation time')
        ax2 = ax1.twinx()
        ax1.set_yscale('log')
        ax1.set_ylabel('permeability')
        ax2.set_ylabel('channelization')
        ax1.plot(t, data[:, 1], 'k', '-', label = f'$\kappa$')
        for i_data in range(3, n_data):
            ax2.plot(t, data[:, i_data], label = f'slice {i_data - 2} / 4')
        ax2.legend()    
        plt.savefig(self.dirname + '/channelization.png')
        plt.close()

    def check_channelization(self, graph: Graph, inc: Incidence, edges: Edges, \
        slice_x: float) -> tuple[int, float]:
        """ Calculate channelization parameter for a slice of the network.

        This function calculates the channelization parameter for a slice of
        the network perpendicular to the main direction of the flow. It checks
        how many edges take half of the total flow going through the slice. The
        function returns the exact number of edges and that number divided by
        the total number of edges in a given slice (so the percentage of edges
        taking half of the total flow in the slice).

        Parameters
        -------
        graph : Graph class object
            network and all its properties

        inc : Incidence class object
            matrices of incidence

        edges : Edges class object
            all edges in network and their parameters

        slice_x : float
            position of the slice

        Returns
        -------
        int
            number of edges taking half of the flow in the slice

        float
            percentage of edges taking half of the flow in the slice
        """
        pos_x = np.array(list(nx.get_node_attributes(graph, 'x').values()))
        # find edges crossing the given slice and their orientation - if edge
        # crosses the slice from left to right, it is marked with 1, if from
        # right to left - -1, if it doesn't cross - 0
        slice_edges = (spr.diags(edges.flow) @ inc.incidence > 0) \
            @ (pos_x <= slice_x) * np.abs(inc.incidence @ (pos_x > slice_x)) \
            - (spr.diags(edges.flow) @ inc.incidence > 0) @ (pos_x > slice_x) \
            * np.abs(inc.incidence @ (pos_x <= slice_x))
        # sort edges from maximum flow to minimum (taking into account
        # their orientation)
        slice_flow = np.array(sorted(slice_edges * edges.fracture_lens \
            * np.abs(edges.flow), reverse = True))
        fraction_flow = 0
        total_flow = np.sum(slice_flow)
        # calculate how many edges take half of the flow
        for i, edge_flow in enumerate(slice_flow):
            fraction_flow += edge_flow
            if fraction_flow > total_flow / 2:
                return (i + 1, (i + 1) / np.sum(slice_flow != 0))
        # if calculation failed, raise an error (it never should happen...)
        raise ValueError("Impossible")
