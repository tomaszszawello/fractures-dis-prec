import data as Da
import dissolution as Di
import draw_net as Dr
import growth as Gr
import incidence as In
import pressure as Pr
import save as Sv

import numpy as np
from build import build
from utils import initialize_iterators, update_iterators
from utils_vtk import save_VTK, save_VTK_nodes
from config import simInputData

sid = simInputData()
G, G2, in_nodes, out_nodes = build(sid)
pressure_b = Pr.create_vector(sid, in_nodes)
cb_b = Di.create_vector(sid, in_nodes)

inc_matrix, mid_matrix, bound_matrix, in_matrix, apertures, fracture_lens, lens, in_edges, \
    out_edges, edge_list = In.create_matrices(sid, G, in_nodes, out_nodes)
iters, tmax, i, t, dt, breakthrough = initialize_iterators(sid)

while t < tmax and i < iters and not breakthrough:
    print(f'Iter {i + 1}/{iters} Time {t:.2f}/{tmax:.2f}')
    pressure, flow = Pr.find_flow(sid, apertures, fracture_lens, lens, inc_matrix, mid_matrix, bound_matrix, in_matrix, pressure_b, in_edges)
    cb = Di.find_cb(sid, apertures, fracture_lens, lens, flow, inc_matrix, in_nodes, out_nodes, cb_b)
    if i % sid.plot_every == 0:
        Da.check_flow(fracture_lens, flow, in_edges, out_edges)
        save_VTK(sid, G, apertures, fracture_lens, lens, np.abs(flow), pressure, cb, name=f'network_{sid.old_iters:04d}.vtk')
        Sv.dump_json_graph(sid, G2, name=f'network_{sid.old_t:04f}')
        G = Pr.update_network(G, edge_list, apertures, np.abs(flow))
        G2 = Pr.update_initial_network(G2, edge_list, sid.b0 * apertures)

    apertures, dt, breakthrough = Gr.update_apertures(sid, flow, cb, apertures, lens, inc_matrix, out_edges, dt)

    i, t = update_iterators(sid, i, t, dt)
    Da.collect_data(sid, pressure)

if i != 1:
    G = Pr.update_network(G, edge_list, apertures, np.abs(flow))
    G2 = Pr.update_initial_network(G2, edge_list, sid.b0 * apertures)
    Da.check_flow(fracture_lens, flow, in_edges, out_edges)
    Sv.dump_json_graph(sid, G2, name=f'network_{sid.old_t:04f}')
    save_VTK(sid, G, apertures, fracture_lens, lens, np.abs(flow), pressure, cb, name=f'network_{sid.old_iters:04d}.vtk')
    save_VTK_nodes(sid, G, in_nodes, out_nodes)
    Da.plot_data(sid)
    Dr.uniform_hist(simInputData, G, in_nodes, out_nodes, apertures, flow, cb, name='d.png', data='d')