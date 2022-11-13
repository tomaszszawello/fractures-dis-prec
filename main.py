import dissolution as Di
import incidence as In
import precipitation as Pi
import pressure as Pr
import save as Sv

from build import build
from utils import check_flow, initialize_iterators, update_iterators, update_diameters_pi
from utils_vtk import save_VTK

from config import simInputData

sid = simInputData()
sid, G, in_nodes, out_nodes, boundary_edges = build(sid)

pressure_b = Pr.create_vector(sid, in_nodes)
cb_b = Di.create_vector(sid, in_nodes)
inc_matrix, mid_matrix, bound_matrix, in_matrix, diams, lens, in_edges, \
    out_edges = In.create_matrices(sid, G, in_nodes, out_nodes)

iters, tmax, i, t, breakthrough = initialize_iterators(sid)

while t < tmax and i < iters and not breakthrough:
    print(f'Iter {i + 1}/{iters} Time {t:.2f}/{tmax:.2f}')

    pressure, flow = Pr.find_flow(sid, diams, lens, inc_matrix, mid_matrix, bound_matrix, in_matrix, pressure_b, in_nodes)
    cb = Di.find_cb(sid, diams, lens, flow, inc_matrix, in_nodes, out_nodes, cb_b)
    cc_b = Pi.create_vector(sid, diams, lens, flow, inc_matrix, in_nodes, cb)
    cc = Pi.find_cc(sid, diams, lens, flow, inc_matrix, in_nodes, out_nodes, cc_b)

    if i % sid.plot_every == 0:
        check_flow(flow, in_edges, out_edges)
        save_VTK(sid, G, boundary_edges, diams, lens, flow, pressure, cb, cc, name=f'network_{sid.old_iters:04d}.vtk')

    # diams, dt = update_diameters(sid, flow, cb, diams, lens, inc_matrix)
    diams, dt, breakthrough = update_diameters_pi(sid, flow, cb, cc, diams, lens, inc_matrix, out_edges)

    i, t = update_iterators(sid, i, t, dt)

if i != 1:
    G = Pr.update_network(G, diams, flow)
    check_flow(flow, in_edges, out_edges)
    save_VTK(sid, G, boundary_edges, diams, lens, flow, pressure, cb, cc, name=f'network_{sid.old_iters:04d}.vtk')
    Sv.save('/save.dill', sid, G, in_nodes, out_nodes, boundary_edges)