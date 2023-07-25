""" Initial parameters of the simulation.

This module contains all parameters set before the simulation. Class
SimInputData is used in nearly all functions. Most of the parameters (apart
from VARIOUS section) are set by the user before starting the simulation.
Most notable parameters are: iters/tmax - simulation length, Da_eff, G - 
dissolution parameters.
"""


# check if pressure & dissolution formulas are correct

class SimInputData:
    """ Configuration class for the whole simulation.
    """
    # GENERAL
    iters: int = 10000
    "maximum number of iterations"
    tmax: float = 500.
    "maximum time"
    plot_every: int = 50
    "frequency (in iterations) of plotting the results"
    collect_every: int = 10
    "frequency (in iterations) of collecting data"
    load_name: str = 'carbonate_x1'
    "name of loaded network"

    # DISSOLUTION & PRECIPITATION
    Da_eff: float = 0.1
    "effective Damkohler number"
    G: float = 1.
    "diffusion to reaction ratio"
    Da: float = Da_eff * (1 + G)
    "Damkohler number"

    # INCLUDE
    include_adt: bool = True
    "include adaptive timestep"
    include_breakthrough: bool = False
    "stop simulation when network is dissolved"

    # INITIAL CONDITIONS
    q_in: float = 1.
    "characteristic flow for inlet edge"
    cb_in: float = 1.
    "inlet B concentration"
    b_break: float = 4.
    "minimal aperture of outlet fracture for network to be dissolved"

    # TIME
    dt: float = 0.01
    "initial timestep (if no adaptive timestep, timestep for whole simulation)"
    growth_rate: float = 0.05
    ("maximum percentage growth of an edges (used for finding adaptive \
     timestep)")
    dt_max: float = 5.
    "maximum timestep (for adaptive)"

    # PARTICLE TRACKING
    track_every = 50
    n_part = 10000
    normalize_channeling = True

    # DRAWING - not sure if necessary
    figsize: float = 10.
    "figure size"
    qdrawconst: float = 15.
    "constant for improving flow drawing"
    ddrawconst: float = 0.5
    "constant for improving diameter drawing"

    # VARIOUS
    n: int = 20
    "size of network (updated later)"
    nsq: int = n ** 2
    "number of nodes (updated later)"
    ne: int = 0
    "number of edges (updated later)"
    b0: float = 1.
    "initial mean aperture (updated later)"
    l0: float = 1.
    "initial mean fracture length (updated later)"
    old_iters: int = 0
    "total iterations of simulation"
    old_t: float = 0.
    "total time of simulation"
    dirname: str = f'G{G:.2f}Daeff{Da_eff:.2f}' + '/' + load_name
    "directory of simulation"
