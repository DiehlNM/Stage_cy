import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot(filename):
    plt.ion()
    plt.rc('axes', labelsize=26)
    plt.rc('axes', titlesize=24)
    data = pd.read_csv(filename)

    mz_max = data['mu_max'].values
    iter = data['iter'].values
    T = data['T'].values
    dim = data['Dim'].values
    energy_per_site = data['energy_per_site'].values

    fig, ax1 = plt.subplots()

    color = 'tab:blue'
    ax1.set_xlabel('Temperature')
    ax1.set_ylabel('Max moment', color=color)
    ax1.plot(T, mz_max, marker='o', color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    """ 
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('mz max', color=color)
    ax2.plot(T, mz_max, marker='x', color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    """

    #plt.title('Energy per site vs temperature')
    fig.tight_layout()
    plt.grid(True)
    plt.ioff()
    plt.show()

plot("MF_Data/List/data_animation.csv")
