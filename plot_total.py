import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot(filename):

    plt.ion()

    data= pd.read_csv(filename)

    mz_max = data['mz_max'].values
    #iter = data['iter'].values
    T = data['T'].values

    plt.plot(T, mz_max, marker='o')
    plt.xlabel('Temperature (T)')
    plt.ylabel('Max moment')
    plt.title('Staggered Magnetization vs Temperature at Half-Filling')
    plt.grid(True)
    plt.ioff()
    plt.show()

plot("MF_Data/List/data_4.csv")