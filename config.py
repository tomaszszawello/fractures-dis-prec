import numpy as np

# zmienić wzory na ciśnienie na d^2, przecałkować wzory na koncentracje przy rozpuszczaniu, zastanowić się, jakie to ma znaczenia dla osadzania, rozważyć wizualizację

class simInputData:
    n = 20 # rozmiar siatki
    iters = 5000 # liczba iteracji
    tmax = 100000
    plot_every = 100
    #save_every = 100

    figsize = 10

    Da_eff = 10
    G = 1

    Da = Da_eff * (1 + G)

    K = 0.6
    Gamma = 1.3

    include_cc = True
    include_vol_a = False


    d0 = 0.100861
    #noise = ["uniform", 1, 0.9] #jednorodny rozkład srednic, srednica początkowa, diameter_wiggle_param
    noise = ["gaussian", 1, 0.1] #gaussowski rozkład srednic, mu, sigma
    #noise = ["lognormal", 1, 0.3] #log-normalny rozkład srednic, mu, sigma

    qin = 1 # przepływ na wejściu
    pout = 0  # jednostki? ciśnienie na wyjściu
    l = 1  # początkowa długosć krawędzi

    cb_in = 1 # mol / dm^3
    cc_in = 0
    vol_a_in = 10

    #alpha = 4 # liczba Sherwooda
    #D = 3e-3 # mm^2/s stała dyfuzji HCl
    #k = 1 # mm/s stała reakcji CaC03 + HCl
    #gamma = 60 # mol/dm^3

    dmin = 0.1
    dmax = 1000

    k_it_th = 0.00001

    adaptive_dt = True
    dt = 0.01
    growth_rate = 0.05 # maksymalny procent średnicy o jaki może urosnąć krawędź
    dt_max = 5

    breakthrough = True
    d_break = 4

    qdrawconst = 15
    ddrawconst = 0.5

    load = 4 # 0 - dane z config, 1 - wczytanie danych z ewoluowanej sieci (plik save), 2 - wczytanie templatki (plik template)
    load_name = 'rect100/G1.00Daeff1.00/7'
    vtk_name = 'network_100x100.vtk'

    #geo = "cylindrical"
    #geo = "donut"
    geo = "own"
    #geo = "top"
    #geo = "own"

    periodic = 'none' #none, top, side, all

    #nodes_own = [[35, 35]]
    #nodes_own = [[40, 40], [60, 60], [60, 40], [40, 60], [50, 64], [50, 36], [36, 50], [64, 50]]
    #in_nodes_own, out_nodes_own = np.array([[25, 50]]) / 100 * n, np.array([[75, 50]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[0, 50]]) / 100 * n, np.array([[100, 50]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50]]) / 100 * n #lista pozycji nodów in i out
    #in_nodes_own, out_nodes_own = np.array([[0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60], [0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50], [5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n
    in_nodes_own, out_nodes_own = np.array([[20, 50]]) / 100 * n, np.array([[80, 50], [70, 25], [70, 75]]) / 100 * n

    nsq = n ** 2
    old_iters = 0
    old_t = 0
    dirname = geo + str(n) + '/' + f'G{G:.2f}Daeff{Da_eff:.2f}'