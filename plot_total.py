import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot(filename):
    plt.ion()

    data = pd.read_csv(filename)

    mz_max = data['mz_max'].values
    iter = data['iter'].values
    T = data['T'].values

    fig, ax1 = plt.subplots()

    color = 'tab:blue'
    ax1.set_xlabel('Temperature (T)')
    ax1.set_ylabel('Max moment', color=color)
    ax1.plot(T, mz_max, marker='o', color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    """
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Number of iterations', color=color)
    ax2.plot(T, iter, marker='x', color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    """

    plt.title('Maximum Magnetization vs Temperature at Half-Filling')
    fig.tight_layout()
    plt.grid(True)
    plt.ioff()
    plt.show()

plot("MF_Data/List/data_4.csv")
