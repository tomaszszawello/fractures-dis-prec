import dissolution as Di
import draw_net as Dr
import incidence as In
import pressure as Pr
import save as Sv

from build import build
from utils import solve_equation, collect_data, update_diameters
from utils_vtk import save_VTK

from config import simInputData
import numpy as np

import scipy.sparse as spr

sid = simInputData()
sid, G, edges, in_nodes, out_nodes, boundary_edges = build(sid) #boundary edges?

presult = Pr.create_vector(sid, in_nodes, out_nodes, edges)
cb_result = Di.create_vector(sid, in_nodes)
inc_matrix, mid_matrix, bound_matrix, in_matrix, diams, lens = In.create_matrices(sid, edges, in_nodes, out_nodes)

iters = sid.old_iters + sid.iters
time = sid.old_t + sid.tmax
t = sid.old_t
i = sid.old_iters
breakthrough = False
while t < time and i < iters and not breakthrough:
    print(f'Iter {i + 1}/{iters} Time {t:.2f}/{time:.2f}')

    #pmatrix2 = Pr.update_matrix(sid, edges, in_nodes, out_nodes)
    #pnow = solve_equation(pmatrix, presult)
    p_matrix = inc_matrix.transpose() @ spr.diags(diams ** 4 / lens) @ inc_matrix
    p_matrix = p_matrix.multiply(mid_matrix) + bound_matrix
    pnow = solve_equation(p_matrix, presult)
    q_in = np.abs(np.sum(diams ** 4 / lens * (in_matrix @ pnow)))
    pnow *= sid.qin * 2 * len(in_nodes) / q_in
    flow = diams ** 4 / lens * (inc_matrix @ pnow)
    cb_inc = np.abs(inc_matrix.transpose() @ (spr.diags(flow) @ inc_matrix > 0))
    qc = flow * np.exp(-sid.Da / (1 + sid.G * diams) * diams * lens / np.abs(flow))
    qc_matrix = np.abs(inc_matrix.transpose() @ spr.diags(qc) @ inc_matrix)
    cb_new = cb_inc.multiply(qc_matrix)
    diag = -np.abs(inc_matrix.transpose()) @ np.abs(flow) / 2
    for node in in_nodes:
        diag[node] = 1
    for node in out_nodes:
        diag[node] *= 2
    cb_new.setdiag(diag)

    #np.savetxt('new.txt', cb_new.toarray())


    #cb_matrix = Di.update_matrix(sid, pnow, edges)
    cb_now = solve_equation(cb_new, cb_result)
    dd = np.abs((spr.diags(flow) @ inc_matrix > 0)) @ cb_now * np.abs(flow)  / (sid.Da * lens * diams) * (1 - np.exp(-sid.Da / (1 + sid.G * diams) * diams * lens / np.abs(flow)))
    

    if i % sid.plot_every == 0:
        G = Pr.update_network(G, sid, edges, diams, flow, in_nodes, out_nodes)
        #Dr.drawvessels(sid, G, in_nodes, out_nodes, boundary_edges, name=f'd_vessels{sid.old_iters:04d}.pdf', data='q')
        #Dr.draw(sid, G, in_nodes, out_nodes, boundary_edges, name=f'q{sid.old_iters:04d}.png', data='q')
        #Dr.draw(sid, G, in_nodes, out_nodes, boundary_edges, name=f'd{sid.old_iters:04d}.png', data='d')
        #Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, pnow, cb_now, name=f'{sid.old_iters:04d}t{t:.2f}.png')
        save_VTK(sid, G, boundary_edges, pnow, cb_now, name=f'network_{sid.old_iters:04d}.vtk')

    if sid.data_collection and i % sid.collect_data_every == 0:
        collect_data(sid, edges, cb_now, pnow, i)

    diams += sid.dt * dd

    #edges, dt, breakthrough = update_diameters(sid, edges, pnow, cb_now)

    if i % sid.save_every == 0 and i != 0:
        Sv.save(f'/save{sid.old_iters:04d}.dill', sid, G, edges, in_nodes, out_nodes, boundary_edges)

    i += 1
    t += sid.dt
    sid.old_iters += 1
    sid.old_t += sid.dt

if i != 1:
    Sv.save('/save.dill', sid, G, edges, in_nodes, out_nodes, boundary_edges)
    G = Pr.update_network(G, sid, edges, diams, flow, in_nodes, out_nodes)
    #Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, pnow, cb_now, name=f'{sid.old_iters:04d}t{t:.2f}.png')
