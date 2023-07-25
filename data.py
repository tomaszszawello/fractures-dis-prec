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

from dataclasses import dataclass
from matplotlib import gridspec
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.sparse as spr

from btc_log import create_bins, create_pdf
from config import SimInputData
from incidence import Edges, Incidence
from network import Graph


class Data():
    """ Contains data collected during the simulation.

    Attributes
    -------
    t : list
        

    pressure : list
        pressure difference between inlet and outlet

    cb_out : list
        difference of inflow and outflow of substance B in the system

    cb_out : list
        difference of inflow and outflow of substance C in the system

    delta_b : float
        current difference of inflow and outflow of substance B in the system

    delta_c : float
        current difference of inflow and outflow of substance C in the system
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
    slices: list = []
    "measures of channeling for slices in the whole system"
    slices_div: list = []
    "measures of channeling for slices in the whole system (divided by ne)"
    slice_names: list = []
    "times of measuring slice channeling"
    tracking: list = []
    "list of breakthrough curves for different times"
    reactive_tracking: list = []
    "list of breakthrough curves with reactive killing for different times"
    concentration_tracking: list = []
    "list of breakthrough curves with concentrations for different times"
    tracking_names: list = []
    "times of tracking"
    velocities: list = []
    "list of collected velocities"
    vol_flow: list = []
    "list of collected volumetric flow"

    def __init__(self, sid: SimInputData, graph: Graph):
        self.dirname = sid.dirname
        self.normalize_channeling = sid.normalize_channeling
        pos_x = np.array(list(nx.get_node_attributes(graph, 'x').values()))
        self.slice_x1 = (np.min(pos_x) + np.average(pos_x)) / 2
        self.slice_x2 = np.average(pos_x)
        self.slice_x3 = (np.max(pos_x) + np.average(pos_x)) / 2

    def save_data(self) -> None:
        """ Save data to text file.

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

    def check_data(self, edges: Edges) -> None:
        """ Check the key physical parameters of the simulation.

        This function calculates and checks if basic physical properties of the
        simulation are valied, i.e. if inflow is equal to outflow.

        Parameters
        -------
        edges : Edges class object
            all edges in network and their parameters
            flow - flow in edges
            inlet - edges connected to inlet nodes
            outlet - edges connected to outlet nodes
        """
        Q_in = np.sum(edges.inlet * np.abs(edges.fracture_lens * edges.flow))
        Q_out = np.sum(edges.outlet * np.abs(edges.fracture_lens * edges.flow))
        print('Q_in =', Q_in, 'Q_out =', Q_out)

    def check_channels(self, graph: Graph, inc: Incidence, edges: Edges, \
        slice_x: np.ndarray):
        pos_x = np.array(list(nx.get_node_attributes(graph, 'x').values()))
        slice_edges = (spr.diags(edges.flow) @ inc.incidence > 0) \
            @ (pos_x <= slice_x) * np.abs(inc.incidence @ (pos_x > slice_x)) \
            - (spr.diags(edges.flow) @ inc.incidence > 0) @ (pos_x > slice_x) \
            * np.abs(inc.incidence @ (pos_x <= slice_x))
        slice_flow = np.array(sorted(slice_edges * edges.fracture_lens \
            * np.abs(edges.flow), reverse = True))
        fraction_flow = 0
        total_flow = np.sum(slice_flow)
        for i, edge_flow in enumerate(slice_flow):
            fraction_flow += edge_flow
            if fraction_flow > total_flow / 2:
                return (i + 1, (i + 1) / np.sum(slice_flow != 0))
        raise ValueError("Impossible")


    def collect_data(self, sid: SimInputData, graph: Graph, inc: Incidence, \
        edges: Edges, p: np.ndarray) -> None:
        """ Collect data from different vectors.

        This function extracts information such as permeability, quantity of
        substances flowing out of the system etc. and saves them in the data
        class.

        Parameters
        -------
        sid : SimInputData class object
            all config parameters of the simulation
            old_t - total time of simulation
            dt - current timestep

        inc : Incidence class object
            matrices of incidence
            incidence - connections of all edges with all nodes

        edges : Edges class object
            all edges in network and their parameters
            flow - flow in edges
            inlet - edges connected to inlet nodes
            outlet - edges connected to outlet nodes

        p : numpy ndarray
            vector of current pressure

        cb : numpy ndarray
            vector of current substance B concentration

        cc : numpy ndarray
            vector of current substance C concentration
        """
        self.t.append(sid.old_t)
        self.perm.append(1 / np.max(p))
        self.order.append((np.sum(edges.flow != 0) - np.sum(np.abs(edges.fracture_lens \
            * edges.flow) ** 2) ** 2 / np.sum(np.abs(edges.fracture_lens \
            * edges.flow) ** 4)) / (np.sum(edges.flow != 0) - 1))
        self.channels_1.append(self.check_channels(graph, inc, edges, \
            self.slice_x1)[1])
        self.channels_2.append(self.check_channels(graph, inc, edges, \
            self.slice_x2)[1])
        self.channels_3.append(self.check_channels(graph, inc, edges, \
            self.slice_x3)[1])

    def plot_data(self) -> None:
        """ Plot data from text file.

        This function loads the data from text file params.txt and plots them
        to file params.png.
        """
        f = open(self.dirname + '/params.txt', 'r', encoding = "utf-8")
        data = np.loadtxt(f)
        n_data = data.shape[1]
        t = data[:, 0]
        plt.figure(figsize = (15, 5))
        plt.suptitle('Parameters')
        spec = gridspec.GridSpec(ncols = n_data - 1, nrows = 1)
        for i_data in range(n_data - 1):
            plt.subplot(spec[i_data]).set_title(f'Data {i_data}')
            plt.plot(t, data[:, i_data + 1])
            plt.xlabel('simulation time')
        plt.savefig(self.dirname + '/params.png')
        plt.close()
        
    def plot_channeling(self) -> None:
        ''' Plot data from text file params.txt and save the plot to params.png.

        Parameters
        -------
        sid : SimInputData class object
            all config parameters of the simulation, here we use attributes:
            dirname - directory of current simulation
        '''
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
        plt.savefig(self.dirname + '/channels.png')
        plt.close()

    def slice_all(self, graph: Graph, inc: Incidence, edges: Edges, name: str):
        pos_x = np.array(list(nx.get_node_attributes(graph, 'x').values()))
        slices = np.linspace(np.min(pos_x), np.max(pos_x), 120)[10:-10]
        channels_tab = []
        channels_tab_div = []
        for x in slices:
            res = self.check_channels(graph, inc, edges, x)
            channels_tab.append(res[0])
            channels_tab_div.append(res[1])
        self.slices.append(channels_tab)
        self.slices_div.append(channels_tab_div)
        self.slice_names.append(name)

    def plot_slices(self, graph: Graph):
        with open(self.dirname + '/slices.txt', "ab") as f:
            np.savetxt(f, self.slices)
            np.savetxt(f, self.slices_div)
        pos_x = np.array(list(nx.get_node_attributes(graph, 'x').values()))
        slices = np.linspace(np.min(pos_x), np.max(pos_x), 120)[10:-10]
        if self.normalize_channeling:
            channeling_0  = np.array(self.slices_div[0])
            for i, channeling in enumerate(self.slices_div):
                plt.plot(slices, np.array(channeling) / channeling_0, \
                    label = self.slice_names[i])
        else:
            for i, channeling in enumerate(self.slices):
                plt.plot(slices, channeling, label = self.slice_names[i])
        plt.xlabel('x')
        plt.ylabel('channeling')
        plt.legend(loc='upper right')
        if self.normalize_channeling:
            plt.savefig(self.dirname + '/slice_norm.png')
        else:
            plt.savefig(self.dirname + '/slice.png')
        plt.close()
        for i, channeling in enumerate(self.slices):
            plt.plot(slices, channeling, label = self.slice_names[i])
        plt.xlabel('x')
        plt.ylabel('channeling')
        plt.legend(loc='upper right')
        plt.savefig(self.dirname + '/slice_no_div.png')
        plt.close()

    def collect_vels(self, edges: Edges):
        """ Save data to text file.

        This function saves the collected data to text file params.txt in
        columns. If the simulation is continued from saved parameters, new data
        is appended to that previously collected.
        """
        self.velocities.append(np.abs(edges.flow / edges.apertures))
        self.vol_flow.append(np.abs(edges.flow * edges.fracture_lens))

    def plot_vels(self) -> None:
        ''' Plot data from text file params.txt and save the plot to params.png.

        Parameters
        -------
        sid : SimInputData class object
            all config parameters of the simulation, here we use attributes:
            dirname - directory of current simulation
        '''
        fig, (ax1, ax2) = plt.subplots(1, 2)
        vel_bx1, vel_pdf1 = create_pdf(data[1], 200)
        vel_bx2, vel_pdf2 = create_pdf(data[3], 200)
        flow_bx1, flow_pdf1 = create_pdf(data[0], 200)
        flow_bx2, flow_pdf2 = create_pdf(data[2], 200)
        ax1.set_title('velocity')
        ax1.loglog(vel_bx1, vel_pdf1,'k',alpha=1, linewidth=2)
        ax1.loglog(vel_bx2, vel_pdf2,'r',alpha=1, linewidth=2)
        ax2.set_title('volumetric flow')
        ax2.loglog(flow_bx1, flow_pdf1,'k',alpha=1, linewidth=2)
        ax2.loglog(flow_bx2, flow_pdf2,'r',alpha=1, linewidth=2)        
        plt.savefig(self.dirname + f'/vels.png')
        plt.close()

    def plot_tracking(self) -> None:
        ''' Plot data from text file params.txt and save the plot to params.png.

        Parameters
        -------
        sid : SimInputData class object
            all config parameters of the simulation, here we use attributes:
            dirname - directory of current simulation
        '''
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (20, 10))
        ax1.set_title('BTC')
        ax2.set_title('BTC with concentration decrease')
        ax3.set_title('BTC with killing particles')
        bx1 = create_bins(self.tracking[-1], 200)
        bx3 = create_bins(self.reactive_tracking[-1], 200)
        for i, time in enumerate(self.tracking_names):
            bx11, pdf = create_pdf(self.tracking[i], 200, x = bx1)
            ax1.loglog(bx11, pdf, alpha = 1, linewidth = 2, label = time)
            bx12, pdf = create_pdf(self.tracking[i], 200, x = bx1, \
                weights = self.concentration_tracking[i])
            ax2.loglog(bx12, pdf, alpha = 1, linewidth = 2, label = time)
            bx13, pdf = create_pdf(self.reactive_tracking[i], 200, x = bx3)
            ax3.loglog(bx13, pdf, alpha = 1, linewidth = 2, label = time)
        ax1.legend()
        ax2.legend()
        ax3.legend()
        plt.savefig(self.dirname + f'/track.png')
        plt.close()


    def find_vels(self, sid: SimInputData, edges: Edges):
        self.tracking_names.append(sid.old_t)
        self.velocities.append(np.abs(edges.flow / edges.apertures))
        self.vol_flow.append(np.abs(edges.flow * edges.fracture_lens))

    def plot_vels(self):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (15, 10))
        ax1.set_title('velocity')
        ax2.set_title('volumetric flow')
        bx1 = create_bins(self.velocities[0], 200)
        bx2 = create_bins(self.vol_flow[0], 200)
        #pdf1, bx1 = np.histogram(self.velocities[0], bins=200, density=False)
        #pdf2, bx2 = np.histogram(self.vol_flow[0], bins=200, density=False)
        for i, time in enumerate(self.tracking_names):
            bx12, pdf = create_pdf(self.velocities[i], 200, x = bx1)
            ax1.loglog(bx12, pdf, alpha = 1, linewidth = 2, label = f'{time:.2f}')
            bx22, pdf = create_pdf(self.vol_flow[i], 200, x = bx2)
            ax2.loglog(bx22, pdf, alpha = 1, linewidth = 2, label = f'{time:.2f}')      
        ax1.legend()
        ax2.legend()
        plt.savefig(self.dirname + f'/vels.png')
        plt.close()