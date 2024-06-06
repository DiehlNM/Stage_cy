import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import Affine2D

def plot_magnetic_moments_3d(filename, dim, vmax, vmin):
    
    plt.ion()
    # Chargement des données
    data = pd.read_csv(filename)

    # Extraire les données mx, my, mz
    mx = data['mx'].values
    my = data['my'].values
    mz = data['mz'].values

    # Initialisation des tableaux pour le réseau carré
    Position = np.zeros((dim * dim, 3))

    # Préparer les données pour le tracé
    x, y, z, u, v, w = [], [], [], [], [], []
    intensitiesxy = np.zeros((dim, dim))
    intensitiesyz = np.zeros((dim, dim))
    intensitiesxz = np.zeros((dim, dim))
    k=0
    for i in range(dim):
        for j in range(dim):
            x.append(i)
            y.append(j)
            idx = i * dim + j
            u.append(mx[idx])
            v.append(my[idx])
            w.append(mz[idx])
            intensitiesxy[i, j] = mx[idx]**2 + my[idx]**2  # Calcul de l'intensité
            intensitiesyz[i, j] = my[idx]**2 + mz[idx]**2
            intensitiesxz[i, j] = mx[idx]**2 + mz[idx]**2

            Position[k, 0] = k  # Numéro de l'atome
            Position[k, 1] = i  # Coordonnée i
            Position[k, 2] = j  # Coordonnée j
            k += 1



    # Préparer les données pour le tracé
    x = Position[:, 1]
    y = Position[:, 2]


    # Plan XY - Vecteurs
    plt.figure(figsize=(10, 10))
    plt.quiver(x, y, u, v)
    plt.xlim([-1, dim])
    plt.ylim([-1, dim])
    plt.title("Magnetic Moments in XY Plane (Vectors)")
    plt.xlabel('X')
    plt.ylabel('Y')
    for i in range(len(x)):
        plt.text(x[i], y[i], str(int(Position[i, 0])), color="blue", fontsize=8)
    plt.show()
    plt.pause(0.001)

    # Plan XY - Intensité
    plt.figure(figsize=(10, 10))
    plt.pcolormesh(np.arange(dim), np.arange(dim), intensitiesxy, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax)
    plt.xlim([-1, dim])
    plt.ylim([-1, dim])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Magnetic Moments in XY Plane (Intensity)")
    plt.colorbar(label='Intensity')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()
    plt.pause(0.001)

    # Plan YZ - Vecteurs
    plt.figure(figsize=(10, 10))
    plt.quiver(x, y, v, w)
    plt.xlim([-1, dim])
    plt.ylim([-1, dim])
    plt.title("Magnetic Moments in YZ Plane (Vectors)")
    plt.xlabel('Y')
    plt.ylabel('Z')
    for i in range(len(y)):
        plt.text(x[i], y[i], str(int(Position[i, 0])), color="blue", fontsize=8)
    plt.show()
    plt.pause(0.001)


    # Plan YZ - Intensité
    plt.figure(figsize=(10, 10))
    plt.pcolormesh(np.arange(dim), np.arange(dim), intensitiesyz, shading='auto', cmap='viridis')
    plt.xlim([-1, dim])
    plt.ylim([-1, dim])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Magnetic Moments in YZ Plane (Intensity)")
    plt.colorbar(label='Intensity')
    plt.xlabel('Y')
    plt.ylabel('Z')
    plt.show()
    plt.pause(0.001)

    # Plan XZ - Vecteurs
    plt.figure(figsize=(10, 10))
    plt.quiver(x, y, u, w)
    plt.xlim([-1, dim])
    plt.ylim([-1, dim])
    plt.title("Magnetic Moments in XZ Plane (Vectors)")
    plt.xlabel('X')
    plt.ylabel('Z')
    for i in range(len(x)):
        plt.text(x[i], y[i], str(int(Position[i, 0])), color="blue", fontsize=8)
    plt.show()
    plt.pause(0.001)

    # Plan XZ - Intensité
    plt.figure(figsize=(10, 10))
    plt.pcolormesh(np.arange(dim), np.arange(dim), intensitiesxz, shading='auto', cmap='viridis')
    plt.xlim([-1, dim])
    plt.ylim([-1, dim])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Magnetic Moments in XZ Plane (Intensity)")
    plt.colorbar(label='Intensity')
    plt.xlabel('X')
    plt.ylabel('Z')
    plt.show()
    plt.pause(0.001)

    plt.ioff()
    plt.show()



# Appel de la fonction avec le chemin vers le fichier CSV et la dimension de la grille
plot_magnetic_moments_3d('Fichier_test/magnetic_data_test20_dim_20_density_0.800000.csv', 20, 0.0732015, 0.00115155)  # Assurez-vous que dim correspond à la dimension de votre grille

#for i in range (7):
#    j = i + 4
#    filename = f"Image/magnetic_data_dim_{j}_density_1.000000.csv"
#    plot_magnetic_moments_3d(filename,j)
