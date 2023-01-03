import delaunay as De
import save as Sv

from geometry import set_geometry
from utils import make_dir
from utils_vtk import build_VTK
from my_network import build_my

from config import simInputData


def build(sid:simInputData):

    if sid.load == 0:
        G, boundary_edges = De.Build_delaunay_net(sid.n, periodic = sid.periodic, noise = sid.noise, dmin = sid.dmin, dmax = sid.dmax)
        in_nodes, out_nodes = set_geometry(G, sid.n, geo=sid.geo, in_nodes=sid.in_nodes_own, out_nodes=sid.out_nodes_own)

        make_dir(sid)

        Sv.save('/template.dill', sid, G, in_nodes, out_nodes, boundary_edges)
        Sv.save_config(sid)

    elif sid.load == 1:
        sid, G, in_nodes, out_nodes, boundary_edges = Sv.load(sid.load_name+'/save.dill')
        # we need to save edges data to be able to load; but save diams and lens will break template

    elif sid.load == 2:
        sid2, G, in_nodes, out_nodes, boundary_edges = Sv.load(sid.load_name+'/template.dill')
        sid.n = sid2.n
        sid.dirname = sid2.dirname + '/template'
        make_dir(sid)
        Sv.save_config(sid)
    
    elif sid.load == 3:
        G, in_nodes, out_nodes, boundary_edges = build_VTK(sid)
        sid.dirname += 'vtk'
        make_dir(sid)
        Sv.save_config(sid)

    elif sid.load == 4:
        G, in_nodes, out_nodes, boundary_edges = build_my(sid)
        sid.dirname += 'my'
        make_dir(sid)
        Sv.save_config(sid)

    return sid, G, in_nodes, out_nodes, boundary_edges