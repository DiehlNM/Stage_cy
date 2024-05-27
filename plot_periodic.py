import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_magnetic_moments_periodicity(filename, dim):
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

    # Création d'une grille 9x9 pour montrer la périodicité
    fig, axs = plt.subplots(3, 3, figsize=(18, 18))

    # Tracer le même graphique dans chaque case de la grille
    for i in range(3):
        for j in range(3):
            axs[i, j].quiver(x, y, u, v)
            axs[i, j].set_xlim([0, dim-1])
            axs[i, j].set_ylim([0, dim-1])
            axs[i, j].set_xticks([])
            axs[i, j].set_yticks([])
            fig.colorbar(axs[i, j].pcolormesh(np.arange(dim), np.arange(dim), intensitiesx, shading='auto', cmap='viridis'), ax=axs[i, j])

    # Ajuster l'espacement pour que les images soient collées
    plt.subplots_adjust(wspace=0, hspace=0)

    # Supprimer les marges autour des subplots
    plt.margins(0, 0)
    
    # Ajouter un titre global
    fig.suptitle("Periodic Magnetic Moments (9x9 Grid)", fontsize=16)
    
    # Ajustement de l'affichage
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

# Appel de la fonction avec le chemin vers le fichier CSV et la dimension de la grille
plot_magnetic_moments_periodicity('Image/magnetic_data_dim_10_density_0.830000.csv', 10)  # Assurez-vous que dim correspond à la dimension de votre grille
