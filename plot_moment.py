import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_magnetic_moments_3d(filename, dim):
    # Chargement des données
    data = pd.read_csv(filename)

    # Extraire les données mx, my, mz
    mx = data['mx'].values
    my = data['my'].values
    mz = data['mz'].values

    # Préparation des données pour le tracé 3D
    x, y, z, u, v, w = [], [], [], [], [], []
    for i in range(dim):
        for j in range(dim):
            # Coordonnées de base pour chaque vecteur
            x.append(i)
            y.append(j)
            z.append(0)  # Utilisation d'un plan z=0 pour simplifier la visualisation
            idx = i * dim + j
            # Composantes des vecteurs
            u.append(mx[idx])
            v.append(my[idx])
            w.append(mz[idx])

    # Création de la figure 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Tracé des vecteurs
    ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)  # Ajustez 'length' selon le besoin

    # Configuration des limites et étiquettes des axes
    ax.set_xlim([0, dim-1])
    ax.set_ylim([0, dim-1])
    ax.set_zlim([-1, 1])  # Ajustez selon la distribution des valeurs de mz

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')


    ax.set_title("3D Magnetic Moments Visualization")

    plt.show()

# Appel de la fonction avec le chemin vers le fichier CSV et la dimension de la grille
plot_magnetic_moments_3d('magnetic_data.csv', 4)  # Assurez-vous que Dim correspond à la dimension de votre grille
