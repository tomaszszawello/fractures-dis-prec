import dissolution as Di
import draw_net as Dr
import growth as Gr
import incidence as In
import pressure as Pr
import save as Sv
import tracking as Tr

import numpy as np
from build import build
from utils import initialize_iterators, update_iterators
from utils_vtk import save_VTK, save_VTK_nodes
from config import SimInputData

sid = SimInputData()
graph, graph_real, inc, edges, data = build(sid)
pressure_b = Pr.create_vector(sid, graph)
cb_b = Di.create_vector(sid, graph)


iters, tmax, i, t, breakthrough = initialize_iterators(sid)

while t < tmax and i < iters and not breakthrough:
    print(f'Iter {i + 1}/{iters} Time {t:.2f}/{tmax:.2f}')
    pressure = Pr.solve_flow(sid, inc, graph, edges, pressure_b)
    cb = Di.solve_dissolution(sid, inc, graph, edges, cb_b)
    if i % sid.collect_every == 0:
        data.collect_data(sid, graph, inc, edges, pressure)
    if int(t) % sid.track_every == 0 and int(t - sid.dt) % sid.track_every != 0:
        #Tr.track(sid, graph, inc, edges, data, pressure)
        data.slice_all(graph, inc, edges, name=f'{sid.old_t:.2f}')
        #data.find_vels(sid, edges)
        #Tr.reactive_kill_track(sid, graph, inc, edges, pressure)
    if i % sid.plot_every == 0:
        data.check_data(edges)
        #graph.update_network(edges)
        graph_real.update_initial_network(sid, edges)
        #save_VTK(sid, G, apertures, fracture_lens, lens, np.abs(flow), pressure, cb, name=f'network_{sid.old_iters:04d}.vtk')
        Sv.dump_json_graph(sid, graph_real, name=f'network_{sid.old_t:04f}')        

    breakthrough, dt = Gr.update_apertures(sid, inc, edges, cb)

    i, t = update_iterators(sid, i, t, dt)

if i != 1:
    #graph.update_network(edges)
    graph_real.update_initial_network(sid, edges)
    #Tr.track(sid, graph, inc, edges, data, pressure)
    #Tr.reactive_kill_track(sid, graph, inc, edges, pressure)
    data.check_data(edges)
    data.slice_all(graph, inc, edges, name=f'{sid.old_t:.2f}')
    #data.find_vels(sid, edges)
    Sv.dump_json_graph(sid, graph_real, name=f'network_{sid.old_t:04f}')
    #save_VTK(sid, G, apertures, fracture_lens, lens, np.abs(flow), pressure, cb, name=f'network_{sid.old_iters:04d}.vtk')
    #save_VTK_nodes(sid, G, in_nodes, out_nodes)
    data.collect_data(sid, graph, inc, edges, pressure)
    data.save_data()
    data.plot_data()
    data.plot_slices(graph)
    #data.plot_tracking()
    #data.plot_vels()
    data.plot_channeling()
    #np.savetxt('vels.txt', data.velocities)
    #Dr.uniform_hist(sid, G, in_nodes, out_nodes, apertures, flow, cb, name='d.png', data='d')
