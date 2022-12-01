import numpy as np
import matplotlib.pyplot as plt

from matplotlib import gridspec

from config import simInputData

def save_data(name, data):
    '''Save data to text file.

    Parameters
    -------
    name : string
        file name

    data : list
        data to append to text file
    '''
    is_saved = False
    while not is_saved: # prevents problems with opening text file
        try:
            f = open(name, 'a')
            np.savetxt(f, [data])
            f.close()
            is_saved = True
        except PermissionError:
            pass

def collect_data(sid:simInputData, pressure):
    '''Collect data from different vectors and save them to text file params.txt.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        dirname - directory of current simulation

    pressure : numpy array
        vector of current pressure
    '''
    data = [np.max(pressure), 0]
    save_data(sid.dirname+'/params.txt', data)

def plot_data(sid:simInputData):
    '''Plot data from text file params.txt and save the plot to params.png.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation, here we use attributes:
        dirname - directory of current simulation
    '''
    f = open(sid.dirname+'/params.txt', 'r')
    data = np.loadtxt(f)
    n_data = data.shape[1]
    len_data = data.shape[0]
    
    plt.figure(figsize=(15, 5))
    plt.suptitle('Parameters')
    spec = gridspec.GridSpec(ncols=n_data, nrows=1)

    for i_data in range(n_data):
        plt.subplot(spec[i_data]).set_title(f'Data {i_data}')
        x = np.linspace(0, len_data, len_data)
        plt.plot(x, data[:, i_data])
        plt.xlabel('iterations')
    
    plt.savefig(sid.dirname +'/params.png')
    plt.close()

def check_flow(flow, in_edges, out_edges):
    ''' Check flow conservation in the system (inflow and outflow)

    Parameters
    -------
    flow : numpy array
        vector of flows in the system
    
    in_edges : numpy array
        list of edges with 1 for inlet edges and 0 otherwise

    out_edges : numpy array
        list of edges with 1 for outlet edges and 0 otherwise

    Returns
    -------
    None
    '''
    Q_in = np.sum(in_edges * np.abs(flow))
    Q_out = np.sum(out_edges * np.abs(flow))
    print('Q_in =', Q_in, 'Q_out =', Q_out)

def check_mass(sid, incidence, flow, cb, vols, in_edges, out_edges, dt, delta_b):
    ''' Check mass balance in the system (inflow of solvent - outflow of solvent = dissolved volume)
    '''
    # np.savetxt('tr.txt', np.abs(incidence.T < 0) @ (np.abs(flow) * in_edges))
    # print (np.sum(np.abs(incidence.T < 0) @ (np.abs(flow) * in_edges)))
    # np.savetxt('tr.txt', np.abs(incidence.T > 0) @ (np.abs(flow) * out_edges))
    # print (np.sum(np.abs(incidence.T > 0) @ (np.abs(flow) * out_edges)))
    delta = (np.abs(incidence.T < 0) @ (np.abs(flow) * in_edges) - np.abs(incidence.T > 0) @ (np.abs(flow) * out_edges)) @ cb * dt
    print (delta)
    delta_b += delta
    vol = np.sum(sid.vol_a_in-vols)
    print (delta_b, sid.Da * vol, delta_b - sid.Da * vol)
    return delta_b