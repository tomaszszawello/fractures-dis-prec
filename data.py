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