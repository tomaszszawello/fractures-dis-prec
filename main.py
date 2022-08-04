#import analysis as An
import draw_net as Dr
import oxygen as Ox
import pressure as Pr
import save as Sv
import upstream as Up
import vegf as Ve

from build import build
from utils import solve_equation, simAnalysisData, collect_data, update_diameters

from config import simInputData
import numpy as np



sid = simInputData()
sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges = build(sid)

presult = Pr.create_vector(sid, in_nodes, out_nodes)
if sid.oxygen:
    bloodoxresult = Ox.create_blood_vector(sid, in_nodes_ox)

iters = sid.old_iters + sid.iters
for i in range(sid.old_iters, iters):
    print(f'Iter {i + 1}/{iters}')

    pmatrix = Pr.update_matrix(sid, edges, in_nodes, out_nodes)
    pnow = solve_equation(pmatrix, presult)
    
    
    if sid.oxygen:
        oxmatrix = Ox.update_matrix(sid, oxresult, bloodoxresult, pnow, edges)
        oxnow = solve_equation(oxmatrix, bloodoxresult)
        vresult = Ve.create_vector(sid, oxresult, oxnow)
    else:
        vresult = Ve.create_vector(sid, oxresult)

    vmatrix = Ve.update_matrix(sid, vresult, edges)
    vnow = solve_equation(vmatrix, vresult)

    if sid.signal:
        smatrix, sresult = Up.update_matrix_upstream(sid, vnow, pnow, edges, in_nodes)
        snow_upstream = solve_equation(smatrix, sresult)
        smatrix, sresult = Up.update_matrix_downstream(sid, vnow, pnow, edges, out_nodes)
        snow_downstream = solve_equation(smatrix, sresult)
        snow = snow_upstream+snow_downstream



    if i % sid.plot_every == 0:
        G = Pr.update_network(G, sid, edges, pnow)

#        Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, oxresult, pnow, oxnow, vnow, snow_upstream, snow, name=f'histogram{sid.old_iters // sid.plot_every:04d}.png')
#        Dr.uniform_hist(sid, G, in_nodes, out_nodes, boundary_edges, oxresult, pnow, pnow, pnow, pnow, pnow, name=f'histogram{sid.old_iters // sid.plot_every:04d}.png')
        #Dr.draw(sid, G, in_nodes, out_nodes, boundary_edges, oxresult, name=f'd{sid.old_iters // sid.save_every:04d}.png', data='d')
#        Dr.drawq(name = f'q{(i+old_iters)//save_every:04d}.png', oxdraw = [])
#        Dr.drawq(name=f'veq{i // save_every:04d}.png', oxdraw=vnow2)
#        Dr.drawq(name=f'oxq{i // save_every:04d}.png', oxdraw=snow)
#        Dr.drawblood(sid, G, in_nodes, out_nodes, boundary_edges, name=f'q_blood{sid.old_iters // sid.save_every:04d}.png', oxresult=oxresult, oxdraw = vnow, data='q')
#        Dr.drawvessels(sid, G, in_nodes, out_nodes, boundary_edges, name=f'q_vessels{sid.old_iters // sid.plot_every:04d}.png', oxresult=oxresult, oxdraw = oxnow, data='q')
        #Dr.drawvessels(sid, G, in_nodes, out_nodes, boundary_edges, name=f'q_vessels{sid.old_iters // sid.plot_every:04d}.png', oxresult=oxresult, oxdraw = oxnow, data='q')
        Dr.drawvessels(sid, G, in_nodes, out_nodes, boundary_edges, name=f'd_vessels{sid.old_iters // sid.plot_every:04d}.pdf', oxresult=oxresult, oxdraw = snow, data='d')
    
    d_pres, d_vegf, d_s = np.zeros(len(edges)), np.zeros(len(edges)), np.zeros(len(edges))
    if sid.shear_d:
        d_pres = Pr.update_graph(sid, edges, pnow)
    if sid.vegf_d:
        d_vegf = Ve.update_graph(sid, vnow, oxresult, edges)
    if sid.signal_d:
        d_s = Up.update_graph(sid, snow_upstream, snow_downstream, pnow, oxresult, edges)
    
    edges = update_diameters(sid, edges, d_pres, d_vegf, d_s)


    oxresult = Ve.update_blood(sid, oxresult, edges)

    if sid.data_collection and i % sid.collect_data_every == 0:
        collect_data(sid, edges, in_nodes, out_nodes, pnow, vnow, oxnow, oxresult)

    if i % sid.save_every == 0 and i != 0:
        Sv.save(f'/save{sid.old_iters}.dill', sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges)

    sid.old_iters += 1
    print (np.average(oxnow ** 2))

Sv.save('/save.dill', sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges)
#if sid.data_collection:
#    Dr.plot_params(sid)

#An.getStrahlerHistogram(G, pnow, oxresult, in_nodes, dirname)


#DiG = An.getDiGraph(G, pnow, oxresult, in_nodes)
#DiG = An.strahlerOrder(DiG)
#An.plotStrahlerGraph(DiG, 'deown101deown_drabina/26')