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
incidence, diams, lens = In.create_incidence_matrix(sid, edges)
pres_inc, in_matrix, diag_matrix, vec_in = In.create_in_out_matrix(sid, edges, in_nodes, out_nodes)
#np.savetxt('pres.txt', pres_inc.toarray())
iters = sid.old_iters + sid.iters
time = sid.old_t + sid.tmax
t = sid.old_t
i = sid.old_iters
breakthrough = False
while t < time and i < iters and not breakthrough:
    print(f'Iter {i + 1}/{iters} Time {t:.2f}/{time:.2f}')

    #pmatrix2 = Pr.update_matrix(sid, edges, in_nodes, out_nodes)
    #pnow = solve_equation(pmatrix, presult)
    # m1 = incidence.transpose() @ diams ** 4 @ lens @ incidence
    # u, s, vh = spr.linalg.svds(m1)
    # pinv = vh.T @ spr.linalg.inv(spr.diags(s)) @ u.T
    # pnow = pres_inc.transpose() @ diams ** 4 @ lens @ pres_inc @ pinv @ presult
    pmatrix = pres_inc.transpose() @ diams ** 4 @ lens @ incidence
    #print (np.sum(diams ** 4 @ lens @ vec_in))
    pmatrix3 = pmatrix.multiply(in_matrix) + np.sum(diams ** 4 @ lens @ vec_in) * diag_matrix
    #np.savetxt('old.txt', pmatrix2.toarray())
    #np.savetxt('new.txt', pmatrix3.toarray())
    #np.savetxt('diag.txt', diag_matrix.toarray())
    #print (type(presult))
    pnow = solve_equation(pmatrix3, presult)
    #np.savetxt('pnow.txt', presult.toarray())


    cb_matrix = Di.update_matrix(sid, pnow, edges)
    cb_now = solve_equation(cb_matrix, cb_result)

    if i % sid.plot_every == 0:
        G = Pr.update_network(G, sid, edges, pnow)
        #Dr.drawvessels(sid, G, in_nodes, out_nodes, boundary_edges, name=f'd_vessels{sid.old_iters:04d}.pdf', data='q')
        #Dr.draw(sid, G, in_nodes, out_nodes, boundary_edges, name=f'q{sid.old_iters:04d}.png', data='q')
        #Dr.draw(sid, G, in_nodes, out_nodes, boundary_edges, name=f'd{sid.old_iters:04d}.png', data='d')
        #Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, pnow, cb_now, name=f'{sid.old_iters:04d}t{t:.2f}.png')
        save_VTK(sid, G, boundary_edges, pnow, cb_now, name=f'network_{sid.old_iters:04d}.vtk')

    if sid.data_collection and i % sid.collect_data_every == 0:
        collect_data(sid, edges, cb_now, pnow, i)

    edges, dt, breakthrough = update_diameters(sid, edges, pnow, cb_now)

    if i % sid.save_every == 0 and i != 0:
        Sv.save(f'/save{sid.old_iters:04d}.dill', sid, G, edges, in_nodes, out_nodes, boundary_edges)

    i += 1
    t += dt
    sid.old_iters += 1
    sid.old_t += dt

if i != 1:
    Sv.save('/save.dill', sid, G, edges, in_nodes, out_nodes, boundary_edges)
    G = Pr.update_network(G, sid, edges, pnow)
    Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, pnow, cb_now, name=f'{sid.old_iters:04d}t{t:.2f}.png')
