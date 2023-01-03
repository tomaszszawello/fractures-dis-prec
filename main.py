import data as Da
import dissolution as Di
import draw_net as Dr
import growth as Gr
import incidence as In
import precipitation as Pi
import pressure as Pr
import save as Sv
import volumes as Vo

from build import build
from utils import initialize_iterators, update_iterators
from utils_vtk import save_VTK

from config import simInputData
import numpy as np

sid = simInputData()
sid, G, in_nodes, out_nodes, boundary_edges = build(sid)
#triangles_inc, triangles_pos = Vo.find_triangles(G)
#vol_a = Vo.create_vector(sid, triangles_pos)
pressure_b = Pr.create_vector(sid, in_nodes)
cb_b = Di.create_vector(sid, in_nodes)
inc_matrix, mid_matrix, bound_matrix, in_matrix, diams, fracture_lens, lens, in_edges, \
    out_edges = In.create_matrices(sid, G, in_nodes, out_nodes)

#delta_b = 0

diams0 = diams.copy()
iters, tmax, i, t, dt, breakthrough = initialize_iterators(sid)

while t < tmax and i < iters and not breakthrough:
    print(f'Iter {i + 1}/{iters} Time {t:.2f}/{tmax:.2f}')

    pressure, flow = Pr.find_flow(sid, diams, fracture_lens, lens, inc_matrix, mid_matrix, bound_matrix, in_matrix, pressure_b, in_nodes)
    cb = Di.find_cb(sid, diams, fracture_lens, lens, flow, inc_matrix, in_nodes, out_nodes, cb_b, dt)
    #cc_b = Pi.create_vector(sid, diams, lens, flow, inc_matrix, in_nodes, cb)
    #cc = Pi.find_cc(sid, diams, lens, flow, inc_matrix, in_nodes, out_nodes, cc_b)
    if i % sid.plot_every == 0:
        Da.check_flow(flow, in_edges, out_edges)
        #save_VTK(sid, G, boundary_edges, diams, lens, flow, pressure, cb, cc, name=f'network_{sid.old_iters:04d}.vtk')
        G = Pr.update_network(G, diams, diams0, flow)
        #Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, cb, cc, vol_a, triangles_pos, name=f'network_{sid.old_iters:.2f}.png')
        Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, cb, name=f'd_{sid.old_iters:.2f}.png', data = 'd')
        Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, cb, name=f'q_{sid.old_iters:.2f}.png', data = 'q')

    diams, vol_a, dt, breakthrough = Gr.update_diameters(sid, flow, cb, diams, fracture_lens, lens, inc_matrix, out_edges, dt)
    #diams, vol_a, dt, breakthrough = Gr.update_diameters(sid, flow, cb, cc, diams, lens, inc_matrix, triangles_inc, vol_a, out_edges, dt)

    i, t = update_iterators(sid, i, t, dt)
    Da.collect_data(sid, pressure)
    #if sid.include_vol_a:
    #    delta_b = Da.check_mass(sid, inc_matrix, flow, cb, vol_a, in_edges, out_edges, dt, delta_b)

if i != 1:
    G = Pr.update_network(G, diams, diams0, flow)
    #Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, cb, cc, vol_a, triangles_pos, name=f'network_{sid.old_iters:.2f}.png')
    Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, cb, name=f'd_{sid.old_iters:.2f}.png', data = 'd')
    Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, cb, name=f'q_{sid.old_iters:.2f}.png', data = 'q')
    Da.check_flow(flow, in_edges, out_edges)
    #save_VTK(sid, G, boundary_edges, diams, lens, flow, pressure, cb, cc, name=f'network_{sid.old_iters:04d}.vtk')
    Sv.save('/save.dill', sid, G, in_nodes, out_nodes, boundary_edges)
    Da.plot_data(sid)