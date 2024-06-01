import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def rotate_vector(u, v, angle):
    """Rotate vectors u, v by a given angle in degrees."""
    theta = np.radians(angle)
    u_rot = u * np.cos(theta) - v * np.sin(theta)
    v_rot = u * np.sin(theta) + v * np.cos(theta)
    return u_rot, v_rot

def plot_magnetic_moments_on_square_lattice(filename, dim):
    # Chargement des données
    data = pd.read_csv(filename)

    # Extraire les données mx, my, mz
    mx = data['mx'].values
    my = data['my'].values
    mz = data['mz'].values

    # Initialisation des tableaux pour le réseau carré
    Position = np.zeros((dim * dim, 3))

    # Remplissage du tableau des positions
    k = 0
    for i in range(dim):
        for j in range(dim):
            Position[k, 0] = k  # Numéro de l'atome
            Position[k, 1] = i  # Coordonnée i
            Position[k, 2] = j  # Coordonnée j
            k += 1

    # Préparer les données pour le tracé
    x = Position[:, 1]
    y = Position[:, 2]
    rotation_angle = 0
    u, v = rotate_vector(mx, my, rotation_angle)  # Appliquer la rotation aux vecteurs dans le plan XY
    u = mx**2
    v = my**2
    w = mz**2

    # Plot XY Plane - Vectors
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

    # Plot XY Plane - Intensity
    plt.figure(figsize=(10, 10))
    intensitiesx = u.reshape((dim, dim))
    plt.pcolormesh(np.arange(dim), np.arange(dim), intensitiesx, shading='auto', cmap='viridis')
    plt.xlim([-1, dim])
    plt.ylim([-1, dim])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Magnetic Moments in XY Plane (Intensity)")
    plt.colorbar(label='Intensity')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

    # Plot YZ Plane - Vectors
    plt.figure(figsize=(10, 10))
    u_yz, w_yz = rotate_vector(my, mz, rotation_angle)  # Appliquer la rotation aux vecteurs dans le plan YZ
    plt.quiver(x, y, u_yz, w_yz)
    plt.xlim([-1, dim])
    plt.ylim([-1, dim])
    plt.title("Magnetic Moments in YZ Plane (Vectors)")
    plt.xlabel('Y')
    plt.ylabel('Z')
    for i in range(len(y)):
        plt.text(x[i], y[i], str(int(Position[i, 0])), color="blue", fontsize=8)
    plt.show()

    # Plot YZ Plane - Intensity
    plt.figure(figsize=(10, 10))
    intensitiesy = v.reshape((dim, dim))
    plt.pcolormesh(np.arange(dim), np.arange(dim), intensitiesy, shading='auto', cmap='viridis')
    plt.xlim([-1, dim])
    plt.ylim([-1, dim])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Magnetic Moments in YZ Plane (Intensity)")
    plt.colorbar(label='Intensity')
    plt.xlabel('Y')
    plt.ylabel('Z')
    plt.show()

    # Plot XZ Plane - Vectors
    plt.figure(figsize=(10, 10))
    u_xz, w_xz = rotate_vector(mx, mz, rotation_angle)  # Appliquer la rotation aux vecteurs dans le plan XZ
    plt.quiver(x, y, u_xz, w_xz)
    plt.xlim([-1, dim])
    plt.ylim([-1, dim])
    plt.title("Magnetic Moments in XZ Plane (Vectors)")
    plt.xlabel('X')
    plt.ylabel('Z')
    for i in range(len(x)):
        plt.text(x[i], y[i], str(int(Position[i, 0])), color="blue", fontsize=8)
    plt.show()

    # Plot XZ Plane - Intensity
    plt.figure(figsize=(10, 10))
    intensitiesz = w.reshape((dim, dim))
    plt.pcolormesh(np.arange(dim), np.arange(dim), intensitiesz, shading='auto', cmap='viridis')
    plt.xlim([-1, dim])
    plt.ylim([-1, dim])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Magnetic Moments in XZ Plane (Intensity)")
    plt.colorbar(label='Intensity')
    plt.xlabel('X')
    plt.ylabel('Z')
    plt.show()

# Appel de la fonction avec le chemin vers le fichier CSV et la dimension de la grille
plot_magnetic_moments_on_square_lattice('Image/magnetic_data10_dim_20_density_1.000000.csv', 20)
