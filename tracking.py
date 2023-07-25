import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as spr

from btc_log import create_pdf
from config import SimInputData
from data import Data
from network import Graph
from incidence import Edges, Incidence

# def track(sid: SimInputData, graph: Graph, inc: Incidence, edges: Edges, \
#     pressure: np.ndarray):
#     print ('Tracking...')
#     breakthrough_times = []
#     neigh_inc = (spr.diags(edges.flow) @ inc.incidence > 0).T
#     tot_flow = np.abs(edges.flow * edges.fracture_lens)
#     tot_time = np.abs(edges.lens / tot_flow)
#     for _ in range(sid.n_part):
#         time = 0
#         inlet_flow = edges.fracture_lens * edges.apertures ** 3 \
#         / edges.lens * (inc.inlet @ p ressure)
#         inlet_flow /= np.sum(inlet_flow)
#         in_edge = np.random.choice(len(inlet_flow), p = inlet_flow)
#         n1, n2 = edges.edge_list[in_edge]
#         if pressure[n1] > pressure[n2]:
#             node = n2
#         else:
#             node = n1
#         while node not in graph.out_nodes:
#             prob = []
#             # get edge indices to neighbors of nodes
#             neigh_edges = neigh_inc[node].nonzero()[1]
#             for edge in neigh_edges:
#                 prob.append(tot_flow[edge])
#             # for neigh in graph.neighbors(node):
#             #     edge_neigh = graph[node][neigh]
#             #     prob.append(np.clip((pressure[node] - pressure[neigh]) * edge_neigh['b'] ** 2, 0, 1000))
#             prob = np.array(prob) / np.sum(prob)
#             edge = neigh_edges[np.random.choice(len(prob), p = prob)]
#             time += tot_time[edge]
#             n1, n2 = edges.edge_list[edge]
#             if n1 == node:
#                 node = n2
#             else:
#                 node = n1
#         breakthrough_times.append(time)
#     plt.hist(breakthrough_times, range = [0, 1000], bins = 50, histtype="step")
#     plt.savefig(sid.dirname + f'/breakthrough_{sid.old_t}.png')
#     plt.close()


def track(sid: SimInputData, graph: Graph, inc: Incidence, \
        edges: Edges, data: Data, pressure: np.ndarray):
    print ('Reactive tracking...')
    breakthrough_times = []
    concentrations = []
    neigh_inc = (spr.diags(edges.flow) @ inc.incidence > 0).T
    tot_flow = np.abs(edges.flow * edges.fracture_lens)
    tot_velocity = np.abs(edges.flow / edges.apertures)
    # tot_time = np.abs(edges.lens / tot_flow)
    tot_time = np.abs(edges.lens / tot_velocity)
    exp = np.exp(-np.abs(2 * sid.Da / (1 + sid.G * edges.apertures) \
        * edges.lens / edges.flow))
    inlet_flow = edges.fracture_lens * edges.apertures ** 3 \
        / edges.lens * (inc.inlet @ pressure)
    inlet_flow /= np.sum(inlet_flow)
    for _ in range(sid.n_part):
        time = 0
        conc = 1
        in_edge = np.random.choice(len(inlet_flow), p = inlet_flow)
        n1, n2 = edges.edge_list[in_edge]
        if pressure[n1] > pressure[n2]:
            node = n2
        else:
            node = n1
        while node not in graph.out_nodes:
            prob = []
            # get edge indices to neighbors of nodes
            neigh_edges = neigh_inc[node].nonzero()[1]
            for edge in neigh_edges:
                prob.append(tot_flow[edge])
            # for neigh in graph.neighbors(node):
            #     edge_neigh = graph[node][neigh]
            #     prob.append(np.clip((pressure[node] - pressure[neigh]) * edge_neigh['b'] ** 2, 0, 1000))
            prob = np.array(prob) / np.sum(prob)
            edge = neigh_edges[np.random.choice(len(prob), p = prob)]
            time += tot_time[edge]
            conc *= exp[edge]
            if conc < 1e-40:
                conc = 0
            n1, n2 = edges.edge_list[edge]
            if n1 == node:
                node = n2
            else:
                node = n1
        breakthrough_times.append(time)
        concentrations.append(conc)

    # bx1, pdf1 = create_pdf(breakthrough_times, 200)
    # bx2, pdf2 = create_pdf(breakthrough_times, 200, weights = concentrations)
    # fig,ax = plt.subplots(figsize=(10,8))
    # ax.loglog(bx1,pdf1,'k',alpha=1,linewidth=2)
    # ax.loglog(bx2,pdf2,'r',alpha=1,linewidth=2)
    # ax.set_xlabel('Time (years)',fontsize=16)
    # ax.set_ylabel('PDF',fontsize=16)
    # #ax.set_title(r'$\alpha = 3.5; \kappa = 100.0$',fontsize=16)
    # ax.tick_params(axis='x', labelsize=16)
    # ax.tick_params(axis='y', labelsize=16)
    # #ax.set_xlim(2e-3,7e-2)
    # #ax.set_ylim(1e-3,1e3)
    # ax.grid(True)
    # fig.savefig(sid.dirname + f'/breakthrough_conc{sid.old_t}.png')
    # plt.close()


    # hist, bin_edges = np.histogram(breakthrough_times, bins = 20)
    # hist_conc, bin_edges = np.histogram(breakthrough_times, bins = 20, weights = concentrations)
    # bins = (bin_edges[:-1] + bin_edges[1:]) / 2
    # logbins = np.logspace(np.log10(1e3), np.log10(1e8), 20)
    # fig, ax1 = plt.subplots()
    # ax1.set_xlabel('breakthrough time')
    # ax1.set_xscale('log')
    # ax1.set_yscale('log')
    # ax1.set_ylim(1e0, 1e4)
    # ax2 = ax1.twinx()
    # ax2.set_yscale('log')
    # ax2.set_ylim(1e-40, 1e4)
    # ax1.set_ylabel('counts')
    # ax2.set_ylabel('reactive counts')
    # ax1.plot(logbins, hist, '.')
    # ax2.plot(logbins, hist_conc, 'r.')
    # plt.savefig(sid.dirname + f'/breakthrough_{sid.old_t}.png')
    # plt.close()

    exp = np.exp(-np.abs(2 * sid.Da / (1 + sid.G * edges.apertures) \
        * edges.lens / edges.flow) / 10)
    breakthrough_times2 = []
    for _ in range(sid.n_part):
        time = 0
        conc = 1
        in_edge = np.random.choice(len(inlet_flow), p = inlet_flow)
        n1, n2 = edges.edge_list[in_edge]
        if pressure[n1] > pressure[n2]:
            node = n2
        else:
            node = n1
        flag = True
        while flag:
            prob = []
            # get edge indices to neighbors of nodes
            neigh_edges = neigh_inc[node].nonzero()[1]
            for edge in neigh_edges:
                prob.append(tot_flow[edge])
            prob = np.array(prob) / np.sum(prob)
            edge = neigh_edges[np.random.choice(len(prob), p = prob)]
            time += tot_time[edge]
            conc *= exp[edge]
            if np.random.rand() > exp[edge]:
                break
            n1, n2 = edges.edge_list[edge]
            if n1 == node:
                node = n2
            else:
                node = n1
            if node in graph.out_nodes:
                flag = False
        if not flag:
            breakthrough_times2.append(time)
            
    data.tracking.append(breakthrough_times)
    data.concentration_tracking.append(concentrations)
    data.reactive_tracking.append(breakthrough_times2)
    data.tracking_names.append(sid.old_t)
    data.velocities.append(np.abs(edges.flow / edges.apertures))
    data.vol_flow.append(np.abs(edges.flow * edges.fracture_lens))

    # is_saved = False
    # while not is_saved: # prevents problems with opening text file
    #     try:
    #         file = open(sid.dirname + '/track1.txt', 'a', \
    #             encoding = "utf-8")
    #         np.savetxt(file, np.array([breakthrough_times, concentrations], \
    #             dtype = float))
    #         file.close()
    #         is_saved = True
    #     except PermissionError:
    #         pass
    # is_saved = False
    # while not is_saved: # prevents problems with opening text file
    #     try:
    #         if sid.old_iters == 0:
    #             file = open(sid.dirname + '/track2.txt', 'a', \
    #                 encoding = "utf-8")
    #         else:
    #             file = open(sid.dirname + '/track3.txt', 'a', \
    #                 encoding = "utf-8")
    #         np.savetxt(file, np.array([breakthrough_times2], dtype = float))
    #         file.close()
    #         is_saved = True
    #     except PermissionError:
    #         pass    
    #bx, pdf1 = create_pdf(breakthrough_times, 200)
    #hist, bin_edges = np.histogram(breakthrough_times, bins = 20)
    
    
    # bx2, pdf2 = create_pdf(breakthrough_times2, 200)
    # fig,ax = plt.subplots(figsize=(10,8))
    # ax.loglog(bx1,pdf1,'k',alpha=1,linewidth=2)
    # ax.loglog(bx2,pdf2,'r',alpha=1,linewidth=2)
    # ax.set_xlabel('Time (years)',fontsize=16)
    # ax.set_ylabel('PDF',fontsize=16)
    # #ax.set_title(r'$\alpha = 3.5; \kappa = 100.0$',fontsize=16)
    # ax.tick_params(axis='x', labelsize=16)
    # ax.tick_params(axis='y', labelsize=16)
    # #ax.set_xlim(2e-3,7e-2)
    # #ax.set_ylim(1e-3,1e3)
    # ax.grid(True)
    # fig.savefig(sid.dirname + f'/breakthrough_kill{sid.old_t}.png')
    # plt.close()
    
    
    
    # bins = (bin_edges[:-1] + bin_edges[1:]) / 2
    # logbins = np.logspace(np.log10(1e3), np.log10(1e8), 20)
    # fig, ax1 = plt.subplots()
    # ax1.set_xlabel('breakthrough time')
    # ax1.set_xscale('log')
    # ax1.set_yscale('log')
    # ax1.set_ylim(1e0, 1e4)
    # ax2 = ax1.twinx()
    # ax2.set_yscale('log')
    # ax2.set_ylim(1e0, 1e4)
    # ax1.set_ylabel('counts')
    # ax2.set_ylabel('reactive counts')
    # ax1.plot(logbins, hist, '.')
    # ax2.plot(logbins, hist_conc, 'r.')
    # plt.savefig(sid.dirname + f'/reactive_{sid.old_t}.png')
    # plt.close()
    # plt.hist(breakthrough_times, range = [0, 1000], bins = 50, histtype="step")
    # plt.savefig(sid.dirname + f'/breakthrough_{sid.old_t}.png')
    # plt.close()
    # plt.hist(breakthrough_times, range = [0, 1000], bins = 50, weights = concentrations, histtype="step")
    # plt.savefig(sid.dirname + f'/breakthrough_conc{sid.old_t}.png')
    # plt.close()




def reactive_kill_track(sid: SimInputData, graph: Graph, inc: Incidence, \
        edges: Edges, pressure: np.ndarray):
    print ('Reactive tracking with killing particles...')
    breakthrough_times = []
    concentrations = []
    neigh_inc = (spr.diags(edges.flow) @ inc.incidence > 0).T
    tot_flow = np.abs(edges.flow * edges.fracture_lens)
    tot_velocity = np.abs(edges.flow / edges.apertures)
    tot_time = np.abs(edges.lens / tot_velocity)
    exp = np.exp(-np.abs(2 * sid.Da / (1 + sid.G * edges.apertures) \
        * edges.lens / edges.flow)) / 10
    inlet_flow = edges.fracture_lens * edges.apertures ** 3 \
        / edges.lens * (inc.inlet @ pressure)
    inlet_flow /= np.sum(inlet_flow)
    while len(breakthrough_times) < sid.n_part:
        time = 0
        conc = 1
        in_edge = np.random.choice(len(inlet_flow), p = inlet_flow)
        n1, n2 = edges.edge_list[in_edge]
        if pressure[n1] > pressure[n2]:
            node = n2
        else:
            node = n1
        flag = True
        while flag:
            prob = []
            # get edge indices to neighbors of nodes
            neigh_edges = neigh_inc[node].nonzero()[1]
            for edge in neigh_edges:
                prob.append(tot_flow[edge])
            # for neigh in graph.neighbors(node):
            #     edge_neigh = graph[node][neigh]
            #     prob.append(np.clip((pressure[node] - pressure[neigh]) * edge_neigh['b'] ** 2, 0, 1000))
            prob = np.array(prob) / np.sum(prob)
            edge = neigh_edges[np.random.choice(len(prob), p = prob)]
            time += tot_time[edge]
            conc *= exp[edge]
            if np.random.rand() > 1 - exp[edge]:
                break
            n1, n2 = edges.edge_list[edge]
            if n1 == node:
                node = n2
            else:
                node = n1
            if node in graph.out_nodes:
                flag = False
        if not flag:
            breakthrough_times.append(time)

    # hist, bin_edges = np.histogram(breakthrough_times, bins = 20)
    # hist_conc, bin_edges = np.histogram(breakthrough_times, bins = 20, weights = concentrations)
    # bins = (bin_edges[:-1] + bin_edges[1:]) / 2
    # logbins = np.logspace(np.log10(1e3), np.log10(1e8), 20)
    # fig, ax1 = plt.subplots()
    # ax1.set_xlabel('breakthrough time')
    # ax1.set_xscale('log')
    # ax1.set_yscale('log')
    # ax1.set_ylim(1e0, 1e4)
    # ax2 = ax1.twinx()
    # ax2.set_yscale('log')
    # ax2.set_ylim(1e-40, 1e4)
    # ax1.set_ylabel('counts')
    # ax2.set_ylabel('reactive counts')
    # ax1.plot(logbins, hist, '.')
    # ax2.plot(logbins, hist_conc, 'r.')
    hist, bin_edges = np.histogram(breakthrough_times, bins = 20)
    bins = (bin_edges[:-1] + bin_edges[1:]) / 2
    logbins = np.logspace(np.log10(1e3), np.log10(1e8), 20)    
    plt.plot(logbins, hist, '.')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(sid.dirname + f'/reactive_{sid.old_t}.png')
    plt.close()
