import delaunay as De
import save as Sv

from geometry import set_geometry, create_edgelist
from oxygen import create_vector
from utils import make_dir, find_ladder

from config import simInputData


def build(sid:simInputData):

    if sid.load == 0:
        G, boundary_edges = De.Build_delaunay_net(sid.n, periodic = sid.periodic, noise = sid.noise, dmin = sid.dmin, dmax = sid.dmax)
        G, in_nodes, out_nodes, reg_nodes, boundary_nodes_out, boundary_nodes_in, in_nodes_ox, out_nodes_ox, oxresult = set_geometry(sid.n, G, geo=sid.geo, R=sid.n//2.5, R_s=sid.n//20, in_nodes=sid.in_nodes_own, out_nodes=sid.out_nodes_own)
        #G, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, oxresult = find_ladder(sid, G)
        edges = create_edgelist(G, in_nodes, out_nodes, reg_nodes, boundary_nodes_out, boundary_nodes_in)

        make_dir(sid)

        Sv.save('/template.dill', sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges)
        Sv.save_config(sid)

    elif sid.load == 1:
        sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges = Sv.load(sid.load_name+'/save.dill')

    elif sid.load == 2:
        sid2, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges = Sv.load(sid.load_name+'/template.dill')
        sid.n = sid2.n
        #sid.dirname = 'params/k/k' + str(sid.k) + 'D' + str(sid.D) + 'Dv' + str(sid.Dv) #sid2.dirname + '/template'
        #sid.dirname = 'params/p' + str(sid.c_pres) + 'v' + str(sid.c_vegf) + 's' + str(sid.c_s)
        #sid.dirname = 'params201/own101/1/template/very-wrong-ox/2'
        sid.dirname = sid2.dirname + '/template'
        make_dir(sid)
        Sv.save_config(sid)

    return sid, G, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges