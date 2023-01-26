import delaunay as De
import save as Sv

from utils import make_dir
from my_network import load_my

from config import simInputData


def build(sid:simInputData):
    G, G2, in_nodes, out_nodes = load_my(sid)
    sid.dirname += 'my'
    make_dir(sid)
    Sv.save_config(sid)

    return G, G2, in_nodes, out_nodes