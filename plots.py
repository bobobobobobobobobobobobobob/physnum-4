import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import os
import subprocess


repertoire = "./"
executable = "physnum_ex4"
input_filename = "configuration.in"


# --- Paramètres physiques ---
epsilon0 = 8.854187817e-12  # permittivité du vide (F/m)
rho0 = epsilon0             # densité volumique de charge (C/m^3)
V0 = 0.0                    # potentiel à r = R (V)
R = 0.05                    # rayon maximal (m)

param = np.linspace(5, 100)
paramstr = 'N1'

N2_coef = 2

fs=16
ls=14



nsimul = param.size

D_data = [None for _ in range(nsimul)]
E_data = [None for _ in range(nsimul)]
phi_data = [None for _ in range(nsimul)]

def _simulate(idx):
    global D_data, E_data, phi_data
    output_file = f"output/{paramstr}={param[idx]}"
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[idx]:.15g} N2={param[idx]*N2_coef:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    D = np.loadtxt(f"{output_file}_D.out")
    E = np.loadtxt(f"{output_file}_E.out")
    phi = np.loadtxt(f"{output_file}_phi.out")
    print(f'Done. ({paramstr}={param[idx]})')
    return D, E, phi

def plot():
    global D_data, E_data, phi_data

    r = np.linspace(0, R, 500)
    # --- Solutions analytiques ---
    phi_ana = V0 + (rho0 * (R**2 - r**2)) / (4 * epsilon0)
    E_ana = (rho0 * r) / (2 * epsilon0)
    D_ana = (rho0 * r) / 2

    r_phi, phi_num = phi_data[0][:, 0], phi_data[0][:, 1]
    r_E, E_num = E_data[0][:, 0], E_data[0][:, 1]
    r_D, D_num = D_data[0][:, 0], D_data[0][:, 1]

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


if __name__ == '__main__':
    batch_size = 10
    for nbatch in range(0, nsimul, batch_size):
        with Pool(processes=batch_size) as pool:
            results = pool.map(_simulate, range(nbatch, min(nsimul, nbatch + batch_size)))
            for i, res in enumerate(results):
                D_data[nbatch + i] = res[0]
                E_data[nbatch + i] = res[1]
                phi_data[nbatch + i] = res[2]
    plot()


