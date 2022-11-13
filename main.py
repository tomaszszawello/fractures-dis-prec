import dissolution as Di
import incidence as In
import pressure as Pr
import save as Sv

from build import build
from utils import initialize_iterators, update_iterators, update_diameters
from utils_vtk import save_VTK

from config import simInputData

sid = simInputData()
sid, G, edges, in_nodes, out_nodes, boundary_edges = build(sid) #boundary edges?

pressure_b = Pr.create_vector(sid, in_nodes)
cb_b = Di.create_vector(sid, in_nodes)
inc_matrix, mid_matrix, bound_matrix, in_matrix, diams, lens = In.create_matrices(sid, edges, in_nodes, out_nodes)
#in_vector, out_vector = In.create_bound_vectors(sid, in_nodes, out_nodes)

iters, tmax, i, t, breakthrough = initialize_iterators(sid)

while t < tmax and i < iters and not breakthrough:
    print(f'Iter {i + 1}/{iters} Time {t:.2f}/{tmax:.2f}')

    pressure, flow = Pr.find_flow(sid, diams, lens, inc_matrix, mid_matrix, bound_matrix, in_matrix, pressure_b, in_nodes)
    cb = Di.find_cb(sid, diams, lens, flow, inc_matrix, in_nodes, out_nodes, cb_b)

    if i % sid.plot_every == 0:
        G = Pr.update_network(G, sid, edges, diams, flow, in_nodes, out_nodes)
        save_VTK(sid, G, boundary_edges, pressure, cb, name=f'network_{sid.old_iters:04d}.vtk')

    diams, dt = update_diameters(sid, flow, cb, diams, lens, inc_matrix)

    # ADD BREAKTHROUGH

    i, t = update_iterators(sid, i, t, dt)

if i != 1:
    G = Pr.update_network(G, sid, edges, diams, flow, in_nodes, out_nodes)
    save_VTK(sid, G, boundary_edges, pressure, cb, name=f'network_{sid.old_iters:04d}.vtk')
    Sv.save('/save.dill', sid, G, edges, in_nodes, out_nodes, boundary_edges)