import numpy as np
import matplotlib.pyplot as plt

# --- Paramètres physiques ---
epsilon0 = 8.854187817e-12  # permittivité du vide (F/m)
rho0 = epsilon0             # densité volumique de charge (C/m^3)
V0 = 0.0                    # potentiel à r = R (V)
R = 0.05                    # rayon maximal (m)

# --- Discrétisation pour les solutions analytiques ---
r = np.linspace(0, R, 500)

# --- Solutions analytiques ---
phi_ana = V0 + (rho0 * (R**2 - r**2)) / (4 * epsilon0)
E_ana = (rho0 * r) / (2 * epsilon0)
D_ana = (rho0 * r) / 2

# === Tracé des solutions analytiques ===
plt.figure()
plt.plot(r, phi_ana, label="Φ(r) analytique", color="blue")
plt.xlabel("r (m)")
plt.ylabel("Φ(r) (V)")
plt.title("Potentiel électrique Φ(r) - Analytique")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(r, E_ana, label="E(r) analytique", color="green")
plt.xlabel("r (m)")
plt.ylabel("E(r) (V/m)")
plt.title("Champ électrique E(r) - Analytique")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(r, D_ana, label="D(r) analytique", color="red")
plt.xlabel("r (m)")
plt.ylabel("D(r) (C/m²)")
plt.title("Déplacement électrique D(r) - Analytique")
plt.grid(True)
plt.legend()
plt.tight_layout()

# === Chargement des solutions numériques ===
fichier_base = "output.out"
phi_data = np.loadtxt(fichier_base + "_phi.out")
E_data = np.loadtxt(fichier_base + "_E.out")
D_data = np.loadtxt(fichier_base + "_D.out")

r_phi, phi_num = phi_data[:, 0], phi_data[:, 1]
r_E, E_num = E_data[:, 0], E_data[:, 1]
r_D, D_num = D_data[:, 0], D_data[:, 1]

# === Tracé des solutions numériques ===
plt.figure()
plt.plot(r_phi, phi_num, label="Φ(r) numérique", color='blue')
plt.xlabel("r (m)")
plt.ylabel("Φ(r) (V)")
plt.title("Potentiel électrique Φ(r) - Numérique")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(r_E, E_num, label="E(r) numérique", color='green')
plt.xlabel("r (m)")
plt.ylabel("E(r) (V/m)")
plt.title("Champ électrique E(r) - Numérique")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(r_D, D_num, label="D(r) numérique", color='red')
plt.xlabel("r (m)")
plt.ylabel("D(r) (C/m²)")
plt.title("Déplacement électrique D(r) - Numérique")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.show()


