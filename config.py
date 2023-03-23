import numpy as np

# zmienić wzory na ciśnienie na d^2, przecałkować wzory na koncentracje przy rozpuszczaniu, zastanowić się, jakie to ma znaczenia dla osadzania, rozważyć wizualizację

class simInputData:
    n = 20 # rozmiar siatki
    iters = 50000 # liczba iteracji
    tmax = 100
    plot_every = 20
    #save_every = 100

    figsize = 10

    Da_eff = 1
    G = 1

    Da = Da_eff * (1 + G)


    b0 = 0.100861
    l0 = 1

    qin = 1 # przepływ na wejściu
    pout = 0  # jednostki? ciśnienie na wyjściu
    l = 1  # początkowa długosć krawędzi

    cb_in = 1 # mol / dm^3


    #alpha = 4 # liczba Sherwooda
    #D = 3e-3 # mm^2/s stała dyfuzji HCl
    #k = 1 # mm/s stała reakcji CaC03 + HCl
    #gamma = 60 # mol/dm^3


    adaptive_dt = True
    dt = 0.01
    growth_rate = 0.05 # maksymalny procent średnicy o jaki może urosnąć krawędź
    dt_max = 5

    breakthrough = False
    b_break = 4

    qdrawconst = 15
    ddrawconst = 0.5

    #load = 4 # 0 - dane z config, 1 - wczytanie danych z ewoluowanej sieci (plik save), 2 - wczytanie templatki (plik template)
    load_name = 'carbonate_x1'
    #vtk_name = 'network_100x100.vtk'

    #geo = "cylindrical"
    #geo = "donut"
    #geo = "own"
    #geo = "top"
    #geo = "own"

    #periodic = 'none' #none, top, side, all

    #nodes_own = [[35, 35]]
    #nodes_own = [[40, 40], [60, 60], [60, 40], [40, 60], [50, 64], [50, 36], [36, 50], [64, 50]]
    #in_nodes_own, out_nodes_own = np.array([[25, 50]]) / 100 * n, np.array([[75, 50]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[0, 50]]) / 100 * n, np.array([[100, 50]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50]]) / 100 * n #lista pozycji nodów in i out
    #in_nodes_own, out_nodes_own = np.array([[0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[40, 40], [60, 60], [60, 40], [40, 60], [0, 50], [15, 85], [50, 100], [85, 85], [100, 50], [85, 15], [50, 0], [15, 15]]) / 100 * n, np.array([[50, 64], [50, 36], [36, 50], [64, 50], [5, 70], [30, 95], [70, 95], [95, 70], [95, 30], [70, 5], [30, 5], [5, 30]]) / 100 * n
    #in_nodes_own, out_nodes_own = np.array([[20, 50]]) / 100 * n, np.array([[80, 50], [70, 25], [70, 75]]) / 100 * n

    nsq = n ** 2
    old_iters = 0
    old_t = 0
    dirname = str(n) + '/' + f'G{G:.2f}Daeff{Da_eff:.2f}'