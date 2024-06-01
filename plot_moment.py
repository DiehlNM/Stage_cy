import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import Affine2D

def plot_magnetic_moments_3d(filename, dim):
    # Chargement des données
    data = pd.read_csv(filename)

    # Extraire les données mx, my, mz
    mx = data['mx'].values
    my = data['my'].values
    mz = data['mz'].values

    # Préparer les données pour le tracé
    x, y, z, u, v, w = [], [], [], [], [], []
    intensitiesx = np.zeros((dim, dim))
    intensitiesy = np.zeros((dim, dim))
    intensitiesz = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            x.append(i)
            y.append(j)
            z.append(0)  # Utilisation d'un plan z=0 pour simplifier la visualisation
            idx = i * dim + j
            u.append(mx[idx])
            v.append(my[idx])
            w.append(mz[idx])
            intensitiesx[i, j] = mx[idx]  # Calcul de l'intensité
            intensitiesy[i, j] = my[idx]
            intensitiesz[i, j] = mz[idx]



    # Création des subplots
    fig, axs = plt.subplots(3, 2)

    transorm = Affine2D().rotate_deg(45)

    # Plan XY - Vecteurs
    axs[0, 0].quiver(x, y, u, v)
    axs[0, 0].set_xlim([0, dim-1])
    axs[0, 0].set_ylim([0, dim-1])
    axs[0, 0].set_title("Magnetic Moments in XY Plane (Vectors)")

    # Plan XY - Intensité
    c1 = axs[0, 1].pcolormesh(np.arange(dim), np.arange(dim), intensitiesx, shading='auto', cmap='viridis')
    #c1.set_transform(transorm + axs[0,1].transData)
    axs[0, 1].set_xlim([0, dim-1])
    axs[0, 1].set_ylim([0, dim-1])
    axs[0, 1].set_title("Magnetic Moments in XY Plane (Intensity)")
    fig.colorbar(c1, ax=axs[0, 1], label='Intensity')

    # Plan YZ - Vecteurs
    axs[1, 0].quiver(x, y, v, w)
    axs[1, 0].set_xlim([0, dim-1])
    axs[1, 0].set_ylim([0, dim-1]) # Ajustez selon la distribution des valeurs de mz
    axs[1, 0].set_title("Magnetic Moments in YZ Plane (Vectors)")

    # Plan YZ - Intensité
    c2 = axs[1, 1].pcolormesh(np.arange(dim), np.arange(dim), intensitiesy, shading='auto', cmap='viridis')
    axs[1, 1].set_xlim([0, dim-1])
    axs[1, 1].set_ylim([0, dim-1])
    axs[1, 1].set_title("Magnetic Moments in YZ Plane (Intensity)")
    fig.colorbar(c2, ax=axs[1, 1], label='Intensity')

    # Plan XZ - Vecteurs
    axs[2, 0].quiver(x, y, u, w)
    axs[2, 0].set_xlim([0, dim-1])
    axs[2, 0].set_ylim([0, dim-1])  # Ajustez selon la distribution des valeurs de mz
    axs[2, 0].set_title("Magnetic Moments in XZ Plane (Vectors)")

    # Plan XZ - Intensité
    c3 = axs[2, 1].pcolormesh(np.arange(dim), np.arange(dim), intensitiesz, shading='auto', cmap='viridis')
    axs[2, 1].set_xlim([0, dim-1])
    axs[2, 1].set_ylim([0, dim-1])
    axs[2, 1].set_title("Magnetic Moments in XZ Plane (Intensity)")
    fig.colorbar(c3, ax=axs[2, 1], label='Intensity')

    # Affichage des plots
    plt.tight_layout()
    plt.show()

# Appel de la fonction avec le chemin vers le fichier CSV et la dimension de la grille
plot_magnetic_moments_3d('Image/magnetic_data10_dim_20_density_1.000000.csv', 20)  # Assurez-vous que dim correspond à la dimension de votre grille

#for i in range (7):
#    j = i + 4
#    filename = f"Image/magnetic_data_dim_{j}_density_1.000000.csv"
#    plot_magnetic_moments_3d(filename,j)
