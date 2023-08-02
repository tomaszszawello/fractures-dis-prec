import dissolution as Di
import growth as Gr
import pressure as Pr
import tracking as Tr

from build import build
from config import SimInputData
from utils import initialize_iterators, update_iterators
from utils_vtk import save_vtk, save_vtk_nodes


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
        Tr.track(sid, inc, graph, edges, data, pressure)
        data.slice_all(graph, inc, edges, name=f'{sid.old_t:.2f}')
    if i % sid.plot_every == 0:
        data.check_data(edges)
        #save_vtk(sid, graph, edges, pressure, cb)
        graph_real.dump_json_graph(sid, edges)

    breakthrough, dt = Gr.update_apertures(sid, inc, edges, cb)

    i, t = update_iterators(sid, i, t, dt)

if i != 1:
    Tr.track(sid, inc, graph, edges, data, pressure)
    #Tr.reactive_kill_track(sid, graph, inc, edges, pressure)
    data.check_data(edges)
    data.slice_all(graph, inc, edges, name=f'{sid.old_t:.2f}')
    #data.find_vels(sid, edges)
    graph_real.dump_json_graph(sid, edges)
    #save_vtk(sid, graph, edges, pressure, cb)
    #save_vtk_nodes(sid, graph)
    data.collect_data(sid, graph, inc, edges, pressure)
    data.save_data()
    data.plot_data()
    data.plot_slices(graph)
    #data.plot_tracking()
    #data.plot_vels()
    data.plot_channeling()
