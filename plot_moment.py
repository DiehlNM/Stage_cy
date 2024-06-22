import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import Affine2D

def plot_magnetic_moments(filename, dim, vmax, vmin):
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
    k = 0
    for i in range(dim):
        for j in range(dim):
            x.append(i)  # Coordonnée x
            y.append(j)  # Coordonnée y
            idx = i * dim + j
            u.append(mx[idx])
            v.append(my[idx])
            w.append(mz[idx])
            intensitiesxy[j, i] = mx[idx]**2 + my[idx]**2  # Calcul de l'intensité
            intensitiesyz[j, i] = my[idx]**2 + mz[idx]**2
            intensitiesxz[j, i] = mx[idx]**2 + mz[idx]**2

            Position[k, 0] = idx  # Numéro de l'atome
            Position[k, 1] = i  # Coordonnée i
            Position[k, 2] = j  # Coordonnée j
            k += 1

    # Préparer les données pour le tracé
    x = Position[:, 1]
    y = Position[:, 2]

    # Figure pour le plan XY
    fig1, axs1 = plt.subplots(1, 2, figsize=(15, 7))

    # Plan XY - Vecteurs
    axs1[0].quiver(x, y, u, v)
    axs1[0].set_xlim([-1, dim])
    axs1[0].set_ylim([-1, dim])
    axs1[0].set_title("Magnetic Moments in XY Plane (Vectors)")
    axs1[0].set_xlabel('X')
    axs1[0].set_ylabel('Y')
    axs1[0].set_xticks(np.arange(-1, dim+1, 1))
    axs1[0].set_yticks(np.arange(-1, dim+1, 1))
    axs1[0].grid(True)

    # Plan XY - Intensité
    im = axs1[1].pcolormesh(np.arange(dim), np.arange(dim), intensitiesxy, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax)
    axs1[1].set_xlim([-1, dim])
    axs1[1].set_ylim([-1, dim])
    axs1[1].set_aspect('equal', adjustable='box')
    axs1[1].set_title("Magnetic Moments in XY Plane (Intensity)")
    fig1.colorbar(im, ax=axs1[1], label='Intensity')
    axs1[1].set_xlabel('X')
    axs1[1].set_ylabel('Y')

    plt.tight_layout()
    plt.show()

    # Figure pour le plan YZ
    fig2, axs2 = plt.subplots(1, 2, figsize=(15, 7))

    # Plan YZ - Vecteurs
    axs2[0].quiver(x, y, v, w)
    axs2[0].set_xlim([-1, dim])
    axs2[0].set_ylim([-1, dim])
    axs2[0].set_title("Magnetic Moments in YZ Plane (Vectors)")
    axs2[0].set_xlabel('Y')
    axs2[0].set_ylabel('Z')
    axs2[0].set_xticks(np.arange(-1, dim+1, 1))
    axs2[0].set_yticks(np.arange(-1, dim+1, 1))
    axs2[0].grid(True)

    # Plan YZ - Intensité
    im = axs2[1].pcolormesh(np.arange(dim), np.arange(dim), intensitiesyz, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax)
    axs2[1].set_xlim([-1, dim])
    axs2[1].set_ylim([-1, dim])
    axs2[1].set_aspect('equal', adjustable='box')
    axs2[1].set_title("Magnetic Moments in YZ Plane (Intensity)")
    fig2.colorbar(im, ax=axs2[1], label='Intensity')
    axs2[1].set_xlabel('Y')
    axs2[1].set_ylabel('Z')
    

    plt.tight_layout()
    plt.show()

    # Figure pour le plan XZ
    fig3, axs3 = plt.subplots(1, 2, figsize=(15, 7))

    # Plan XZ - Vecteurs
    axs3[0].quiver(x, y, u, w)
    axs3[0].set_xlim([-1, dim])
    axs3[0].set_ylim([-1, dim])
    axs3[0].set_title("Magnetic Moments in XZ Plane (Vectors)")
    axs3[0].set_xlabel('X')
    axs3[0].set_ylabel('Z')
    axs3[0].set_xticks(np.arange(-1, dim+1, 1))
    axs3[0].set_yticks(np.arange(-1, dim+1, 1))
    axs3[0].grid(True)

    # Plan XZ - Intensité
    im = axs3[1].pcolormesh(np.arange(dim), np.arange(dim), intensitiesxz, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax)
    axs3[1].set_xlim([-1, dim])
    axs3[1].set_ylim([-1, dim])
    axs3[1].set_aspect('equal', adjustable='box')
    axs3[1].set_title("Magnetic Moments in XZ Plane (Intensity)")
    fig3.colorbar(im, ax=axs3[1], label='Intensity')
    axs3[1].set_xlabel('X')
    axs3[1].set_ylabel('Z')

    plt.tight_layout()
    plt.show()

# Appel de la fonction avec le chemin vers le fichier CSV et la dimension de la grille
plot_magnetic_moments('MF_Data/Data_2/magnetic_data_dim_21_density_1.000000_temperature0.000001.csv', 20, 0.08, 0.0)  # Assurez-vous que dim correspond à la dimension de votre grille
