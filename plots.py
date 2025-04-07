import numpy as np
import matplotlib.pyplot as plt

# Spécifie le nom de base du fichier de sortie utilisé dans le fichier C++
# Exemple : si fichier="result", il doit y avoir "result_phi.out", etc.
fichier_base = "output.out"  # À adapter si besoin

# Chargement des données
phi_data = np.loadtxt(fichier_base + "_phi.out")
E_data = np.loadtxt(fichier_base + "_E.out")
D_data = np.loadtxt(fichier_base + "_D.out")

r_phi, phi = phi_data[:, 0], phi_data[:, 1]
r_E, E = E_data[:, 0], E_data[:, 1]
r_D, D = D_data[:, 0], D_data[:, 1]

# Création des figures
plt.figure(figsize=(12, 8))

# Potentiel électrique
plt.subplot(3, 1, 1)
plt.plot(r_phi, phi, label="Potentiel électrique $\\phi(r)$", color='blue')
plt.xlabel("Rayon r (m)")
plt.ylabel("$\\phi$ (V)")
plt.grid(True)
plt.legend()

# Champ électrique
plt.subplot(3, 1, 2)
plt.plot(r_E, E, label="Champ électrique $E(r)$", color='red')
plt.xlabel("Rayon r (m)")
plt.ylabel("$E$ (V/m)")
plt.grid(True)
plt.legend()

# Champ de déplacement électrique
plt.subplot(3, 1, 3)
plt.plot(r_D, D, label="Déplacement électrique $D(r)$", color='green')
plt.xlabel("Rayon r (m)")
plt.ylabel("$D$ (C/m²)")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

