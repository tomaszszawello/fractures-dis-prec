import dissolution as Di
import draw_net as Dr
import pressure as Pr
import save as Sv

from build import build
from utils import solve_equation, collect_data, update_diameters

from config import simInputData
import numpy as np



sid = simInputData()
sid, G, edges, in_nodes, out_nodes, boundary_edges = build(sid) #boundary edges?

presult = Pr.create_vector(sid, in_nodes, out_nodes, edges)
cb_result = Di.create_vector(sid, in_nodes)



iters = sid.old_iters + sid.iters
for i in range(sid.old_iters, iters):
    print(f'Iter {i + 1}/{iters}')

    pmatrix = Pr.update_matrix(sid, edges, in_nodes, out_nodes)
    pnow = solve_equation(pmatrix, presult)

    cb_matrix = Di.update_matrix(sid, pnow, edges)
    cb_now = solve_equation(cb_matrix, cb_result)

    if i % sid.plot_every == 0:
        G = Pr.update_network(G, sid, edges, pnow)
        #Dr.drawvessels(sid, G, in_nodes, out_nodes, boundary_edges, name=f'd_vessels{sid.old_iters:04d}.pdf', data='q')
        #Dr.draw(sid, G, in_nodes, out_nodes, boundary_edges, name=f'q{sid.old_iters:04d}.png', data='q')
        #Dr.draw(sid, G, in_nodes, out_nodes, boundary_edges, name=f'd{sid.old_iters:04d}.png', data='d')
        Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, pnow, cb_now, name=f'{sid.old_iters:04d}.png')


    #d_pres = Pr.update_graph(sid, edges, pnow)
    edges = update_diameters(sid, edges, pnow, cb_now)

    if sid.data_collection and i % sid.collect_data_every == 0:
        collect_data(sid, edges, in_nodes, out_nodes, pnow)

    if i % sid.save_every == 0 and i != 0:
        Sv.save(f'/save{sid.old_iters:04d}.dill', sid, G, edges, in_nodes, out_nodes, boundary_edges)

    dlist = []
    for e in edges:
        (n1, n2, d, l, t) = e
        dlist.append(d)

    sid.old_iters += 1

Sv.save('/save.dill', sid, G, edges, in_nodes, out_nodes, boundary_edges)