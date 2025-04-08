import numpy as np
import matplotlib.pyplot as plt

# --- Paramètres physiques ---
epsilon0 = 8.854187817e-12  # permittivité du vide (F/m)
rho0 = epsilon0             # densité volumique de charge (C/m^3)
V0 = 0.0                    # potentiel à r = R (V)
R = 0.05                    # rayon maximal (m)

fs=16
ls=14

# --- Discrétisation pour les solutions analytiques ---
r = np.linspace(0, R, 500)

# --- Solutions analytiques ---
phi_ana = V0 + (rho0 * (R**2 - r**2)) / (4 * epsilon0)
E_ana = (rho0 * r) / (2 * epsilon0)
D_ana = (rho0 * r) / 2

# === Chargement des solutions numériques ===
fichier_base = "output.out"
phi_data = np.loadtxt(fichier_base + "_phi.out")
E_data = np.loadtxt(fichier_base + "_E.out")
D_data = np.loadtxt(fichier_base + "_D.out")

r_phi, phi_num = phi_data[:, 0], phi_data[:, 1]
r_E, E_num = E_data[:, 0], E_data[:, 1]
r_D, D_num = D_data[:, 0], D_data[:, 1]

# === Tracé des solutions numériques et analytiques ===
plt.figure()
plt.plot(r_phi, phi_num, label="Solution Numérique", color='blue')
plt.plot(r, phi_ana, label="Solution Analytique", color="violet")
plt.xlabel(r"$r$ [m]",fontsize=fs)
plt.ylabel(r"$\Phi$(r) [V]",fontsize=fs)
plt.grid(True)
plt.ticklabel_format(style='scientific', scilimits=(0,0))
plt.legend(fontsize=ls)
plt.tight_layout()

plt.figure()
plt.plot(r_E, E_num, label="Solution Numérique", color='blue')
plt.plot(r, E_ana, label="Solution Analytique", color="violet")
plt.xlabel(r"$r$ [m]",fontsize=fs)
plt.ylabel(r"$E(r)$ [V/m]",fontsize=fs)
plt.grid(True)
plt.ticklabel_format(style='scientific', scilimits=(0,0))
plt.legend(fontsize=ls)
plt.tight_layout()


plt.figure()
plt.plot(r_D, D_num, label="Solution numérique", color='blue')
plt.plot(r, D_ana, label="Solution analytique", color="violet")
plt.xlabel(r"$r$ [m]",fontsize=fs)
plt.ylabel(r"$D(r)$ [C/m²]",fontsize=fs)
plt.grid(True)
plt.ticklabel_format(style='scientific', scilimits=(0,0))
plt.legend(fontsize=ls)
plt.tight_layout()


# --- Discrétisation cylindrique ---
r_vals = np.linspace(0, R, 200)
theta_vals = np.linspace(0, 2*np.pi, 200)
r, theta = np.meshgrid(r_vals, theta_vals)

# --- Coordonnées cartésiennes ---
X = r * np.cos(theta)
Y = r * np.sin(theta)

# --- Potentiel analytique Φ(r) ---
phi = V0 + (rho0 * (R**2 - r) ** 2) / (4 * epsilon0)

# --- Affichage 3D ---
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(X, Y, phi, cmap='plasma', edgecolor='none')
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
#ax.set_zlabel(r"$\Phi$ [V]")
fig.colorbar(surf, ax=ax, label=r"$\Phi$ [V]")
plt.show()

