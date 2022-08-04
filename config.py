import numpy as np
from utils import fParams

class simInputData:
    n = 201 # rozmiar siatki
    iters = 1  # liczba iteracji
    plot_every = 1
    save_every = 200
    collect_data_every =  30

    figsize = 10

    #noise = ["uniform", 1, 0.9] #jednorodny rozkład srednic, srednica początkowa, diameter_wiggle_param
    noise = ["gaussian", 1, 0.2] #gaussowski rozkład srednic, mu, sigma
    #noise = ["lognormal", 1, 0.3] #log-normalny rozkład srednic, mu, sigma

    oxygen = True
    signal = True

    shear_d = True
    vegf_d = True
    signal_d = True

    data_collection = True

    qin = 20 * (n / 201) # ilosć wpływającej krwi
    presout = 0  # cisnienie na wyjsciu
    mu = 0.0035  # współczynnik lepkosci
    l = 1  # początkowa długosć krawędzi
    c1 = np.pi / 128  # stała przepływu
    c2 = 64 / np.pi  # stała siły


    dth = 5 # graniczna grubosć
    dmin = 0.01
    dmax = 20

    Dv = 0.2 * (n / 201) # współczynnik dyfuzji VEGF
    D_ox_b = 2e-5 * (n / 201) # współczynnik dyfuzji tlenu w krwi
    k_ox_b = 1e-3 / qin * (n / 201) # współczynnik reakcji tlenu w krwi
    D_ox_t = 2e-3 * (n / 201) # współczynnik dyfuzji tlenu w tkance
    k_ox_t = 7e-6 * (n / 201) # współczynnik reakcji tlenu w tkance

    linear = 1

    c_pres = 0.1
    c_vegf = 1
    c_s = 0.05

    decrease = 0

    F_p = fParams({'F0': 0.0002, 'F1': 0.002, 'z0': 0, 'z1': 0.03})
    F_ox = fParams({'F0': 20, 'F1': 300, 'z0': 0, 'z1': 0.4})
    F_s = fParams({'F0': 10000, 'F1': 4000000, 'z0': 0, 'z1': 0.5})
    
    R_c = 0.3
    R_a = 15


    ks = 0.05
    v = -1 * (n / 201)
    R = 1

    cox_in = 1

    q_prun = 0.02

    qdrawconst = 15
    ddrawconst = 0.3

    load = 2 # 0- dane z config, 1- wczytanie danych z ewoluowanej sieci (plik save), 2- wczytanie templatki (plik template)
    #load_name = 'own101/117/template/13'
    load_name = 'own201/15'


    #geo = "cylindrical"
    #geo = "donut"
    #geo = "rect"
    #geo = "top"
    geo = "own"

    periodic = 'all' #none, top, side, all

    #nodes_own = [[35, 35]]
    #nodes_own = [[40, 40], [60, 60], [60, 40], [40, 60], [50, 64], [50, 36], [36, 50], [64, 50]]
    in_nodes_own, out_nodes_own = np.array([[25, 50]]) / 100 * n, np.array([[75, 50]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[0, 50]]) / 100 * n, np.array([[100, 50]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50]]) / 100 * n #lista pozycji nodów in i out
    #in_nodes_own, out_nodes_own = np.array([[0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60], [0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50], [5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n
    
    nsq = n ** 2
    old_iters = 0
    dirname = geo + str(n)