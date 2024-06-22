import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from matplotlib.animation import FuncAnimation


def plot_magnetic_moments(axs, filename, dim, vmax, vmin):
    # Chargement des données
    data = pd.read_csv(filename)

    # Extraire les données mx, my, mz
    mx = data['mx'].values
    my = data['my'].values
    mz = data['mz'].values

    # Préparer les données pour le tracé
    x, y, z, u, v, w = [], [], [], [], [], []
    intensitiesxy = np.zeros((dim, dim))
    intensitiesyz = np.zeros((dim, dim))
    intensitiesxz = np.zeros((dim, dim))
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

    # Mise à jour des graphiques pour le plan XY
    axs[0][0].cla()
    axs[0][0].quiver(x, y, u, v)
    axs[0][0].set_xlim([-1, dim])
    axs[0][0].set_ylim([-1, dim])
    axs[0][0].set_title("Magnetic Moments in XY Plane (Vectors)")
    axs[0][0].set_xlabel('X')
    axs[0][0].set_ylabel('Y')
    axs[0][0].grid(True)

    axs[0][1].cla()
    im = axs[0][1].pcolormesh(np.arange(dim), np.arange(dim), intensitiesxy, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax)
    axs[0][1].set_xlim([-1, dim])
    axs[0][1].set_ylim([-1, dim])
    axs[0][1].set_aspect('equal', adjustable='box')
    axs[0][1].set_title("Magnetic Moments in XY Plane (Intensity)")
    axs[0][1].set_xlabel('X')
    axs[0][1].set_ylabel('Y')

    # Mise à jour des graphiques pour le plan YZ
    axs[1][0].cla()
    axs[1][0].quiver(x, y, v, w)
    axs[1][0].set_xlim([-1, dim])
    axs[1][0].set_ylim([-1, dim])
    axs[1][0].set_title("Magnetic Moments in YZ Plane (Vectors)")
    axs[1][0].set_xlabel('Y')
    axs[1][0].set_ylabel('Z')
    axs[1][0].grid(True)

    axs[1][1].cla()
    im = axs[1][1].pcolormesh(np.arange(dim), np.arange(dim), intensitiesyz, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax)
    axs[1][1].set_xlim([-1, dim])
    axs[1][1].set_ylim([-1, dim])
    axs[1][1].set_aspect('equal', adjustable='box')
    axs[1][1].set_title("Magnetic Moments in YZ Plane (Intensity)")
    axs[1][1].set_xlabel('Y')
    axs[1][1].set_ylabel('Z')

    # Mise à jour des graphiques pour le plan XZ
    axs[2][0].cla()
    axs[2][0].quiver(x, y, u, w)
    axs[2][0].set_xlim([-1, dim])
    axs[2][0].set_ylim([-1, dim])
    axs[2][0].set_title("Magnetic Moments in XZ Plane (Vectors)")
    axs[2][0].set_xlabel('X')
    axs[2][0].set_ylabel('Z')
    axs[2][0].grid(True)

    axs[2][1].cla()
    im = axs[2][1].pcolormesh(np.arange(dim), np.arange(dim), intensitiesxz, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax)
    axs[2][1].set_xlim([-1, dim])
    axs[2][1].set_ylim([-1, dim])
    axs[2][1].set_aspect('equal', adjustable='box')
    axs[2][1].set_title("Magnetic Moments in XZ Plane (Intensity)")
    axs[2][1].set_xlabel('X')
    axs[2][1].set_ylabel('Z')

    plt.tight_layout()

# Répertoire contenant les fichiers CSV
directory = "MF_Data/Data_2/"
file_list = glob.glob(os.path.join(directory, "*.csv"))

# Assurez-vous que les fichiers sont triés dans l'ordre voulu
file_list.sort()

# Dimensions et autres paramètres (à ajuster selon vos besoins)
dim = 4
vmax = 1.0
vmin = -1.0

fig, axs = plt.subplots(3, 2, figsize=(15, 21))

def update_plot(i):
    plot_magnetic_moments(axs, file_list[i], dim, vmax, vmin)

# Créer l'animation
ani = FuncAnimation(fig, update_plot, frames=len(file_list), repeat=False)

# Sauvegarder l'animation ou afficher
ani.save('MF_Data/Animation/magnetic_moments_animation.mp4', writer='ffmpeg', fps=2)
plt.show()  # Si vous voulez afficher l'animation

