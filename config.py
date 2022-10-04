import numpy as np
from utils import fParams

class simInputData:
    n = 201 # rozmiar siatki
    iters = 901  # liczba iteracji
    plot_every = 50
    save_every = 1000
    collect_data_every =  1000

    figsize = 10

    #noise = ["uniform", 1, 0.9] #jednorodny rozkład srednic, srednica początkowa, diameter_wiggle_param
    noise = ["gaussian", 0.1, 0.02] #gaussowski rozkład srednic, mu, sigma
    #noise = ["lognormal", 1, 0.3] #log-normalny rozkład srednic, mu, sigma

    shear_d = False

    data_collection = False

    qin = 0.03 # jednostki? przepływ na wejściu
    #pin = 1 # jednostki? ciśnienie na wejściu
    pout = 0  # jednostki? ciśnienie na wyjściu
    mu = 1e-3 #kg/(m s) współczynnik lepkości dynamicznej
    l = 1  # początkowa długosć krawędzi
    c1 = np.pi / 128  # stała przepływu


    cb_in = 4 # mol / dm^3

    alpha = 4 # liczba Sherwooda
    D = 3e-3 # mm^2/s stała dyfuzji HCl
    k = 1 # mm/s stała reakcji CaC03 + HCl
    gamma = 60 # mol/dm^3
    dt = 100 # s krok czasowy



    dmin = 0.01
    dmax = 20


    qdrawconst = 15
    ddrawconst = 1

    load = 0 # 0- dane z config, 1- wczytanie danych z ewoluowanej sieci (plik save), 2- wczytanie templatki (plik template)
    load_name = 'rect101/1'


    #geo = "cylindrical"
    #geo = "donut"
    geo = "rect"
    #geo = "top"
    #geo = "own"

    periodic = 'top' #none, top, side, all

    #nodes_own = [[35, 35]]
    #nodes_own = [[40, 40], [60, 60], [60, 40], [40, 60], [50, 64], [50, 36], [36, 50], [64, 50]]
    in_nodes_own, out_nodes_own = np.array([[25, 50]]) / 100 * n, np.array([[75, 50]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[0, 50]]) / 100 * n, np.array([[100, 50]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50]]) / 100 * n #lista pozycji nodów in i out
    #in_nodes_own, out_nodes_own = np.array([[0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60], [0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50], [5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n
    
    nsq = n ** 2
    old_iters = 0
    G = k * noise[1] / D / alpha
    Da = np.pi * noise[1] * k * l / qin / (1 + G)
    dirname = geo + str(n) + '/' + f'G{G:.2f}Da{Da:.2f}'