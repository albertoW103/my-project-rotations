#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import math


def get_blocks_proteins(input_filename):
    def _to_float(s):
        return float(s.replace('D', 'e').replace('d', 'e'))

    blocks = []
    with open(input_filename, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if line == "":
                continue
            try:
                n = int(line)
            except ValueError:
                continue

            header = f.readline()
            if not header:
                break

            frame = []
            for _ in range(n):
                atom_line = f.readline()
                if not atom_line:
                    raise ValueError("XYZ truncado al leer coordenadas.")
                parts = atom_line.split()
                if len(parts) < 4:
                    raise ValueError(f"Línea inválida en XYZ: {atom_line!r}")
                res = parts[0]
                x = _to_float(parts[1]); y = _to_float(parts[2]); z = _to_float(parts[3])
                frame.append([res, x, y, z])

            blocks.append(frame)
    return blocks


def get_angles(xrot, yrot, zrot, n1, n2, n3):
    """
    Devuelve (cosθ, φ [rad], cos²θ) del vector n1→n2 (índices 1-based).
    """
    if n1 is None or n2 is None:
        return None
    nu = len(xrot)

    xrot = np.asarray(xrot, float)
    yrot = np.asarray(yrot, float)
    zrot = np.asarray(zrot, float)
    
    # ordena los indices de xyz to python
    n1 -= 1; n2 -= 1; n3 -= 1
    
    # comprobacion:
    if not (0 <= n1 < nu and 0 <= n2 < nu and 0 <= n3 < nu):
        raise IndexError("n1/n2/n3 out of range.")

    delx = xrot[n1] - xrot[n2]
    dely = yrot[n1] - yrot[n2]
    delz = zrot[n1] - zrot[n2]
    
    delr = math.sqrt(delx*delx + dely*dely + delz*delz)
    
    costh = delz / delr
    phi   = math.atan2(dely, delx)
    costh2 = costh*costh
    
    # --- ψ: ángulo orientado en [0, 2π) en convención ZYZ, usando solo trigonometría ---

    # Vector v13 = r(n1) - r(n3)
    delx13 = xrot[n1] - xrot[n3]
    dely13 = yrot[n1] - yrot[n3]
    delz13 = zrot[n1] - zrot[n3]
    
    # sen(θ) y cos(θ) a partir de cosθ (ya calculado)
    sinth = math.sqrt(max(0.0, 1.0 - costh*costh))
    cosph = math.cos(phi)
    sinph = math.sin(phi)
    
    # Ejes del cuerpo tras aplicar φ y θ (columnas 1 y 2 de Rz(φ)·Ry(θ)):
    # x'' = (cosφ cosθ,  sinφ cosθ, -sinθ)
    exx = cosph * costh
    exy = sinph * costh
    exz = -sinth
    
    # y'' = (-sinφ, cosφ, 0)
    eyx = -sinph
    eyy =  cosph
    eyz =  0.0
    
    # Proyección de v13 sobre (x'', y''):
    px = delx13*exx + dely13*exy + delz13*exz
    py = delx13*eyx + dely13*eyy + delz13*eyz
    
    # Ángulo orientado en el plano x''–y'':
    psi = math.atan2(py, px)
           
    return costh, phi, costh2, psi


def plot_angles_1x3(costheta_list, phi_list, psi_list, outname):

    outname="histos_angles_1x3.png"
    
    costheta_list = np.array(costheta_list)
    phi_list      = np.array(phi_list) / np.pi   # [0,2]
    psi_list      = np.array(psi_list) / np.pi   # [0,2]
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    bins=80
    
    # cos(theta) en [-1, 1]
    ax = axes[0]
    ax.hist(costheta_list, bins=bins, range=(-1, 1), density=True, alpha=0.8, edgecolor='black')
    ax.axhline(0.5, linestyle='dashed', color='black')
    ax.set_title('Density of $\cos\\theta$')
    ax.set_ylabel('Density')
    ax.set_xlabel(r'$\cos\theta$')
    ax.set_ylim(0, 1)
    ax.set_xlim(-1, 1)
    ax.grid(True, alpha=0.3)
    
    # φ en [0, 2π)
    ax = axes[1]
    ax.hist(phi_list, bins=bins, range=(-1, 1), density=True, alpha=0.8, edgecolor='black')
    ax.axhline(0.5, linestyle='dashed', color='black')
    ax.set_title('Density of $\phi$')
    ax.set_ylabel('Density')
    ax.set_xlabel(r'$\phi$/$\pi$')
    ax.set_ylim(0, 1)
    ax.set_xlim(-1, 1)
    ax.grid(True, alpha=0.3)
    
    # ψ en [0, 2π)
    ax = axes[2]
    ax.hist(psi_list, bins=bins, range=(-1, 1), density=True, alpha=0.8, edgecolor='black')
    ax.axhline(0.5, linestyle='dashed', color='black')
    ax.set_title('Density of $\psi$')
    ax.set_ylabel('Density')
    ax.set_xlabel(r'$\psi$/$\pi$')
    ax.set_ylim(0, 1)
    ax.set_xlim(-1, 1)
    ax.grid(True, alpha=0.3)

    plt.savefig(outname, dpi=300, bbox_inches='tight')
    print(f"Saved: {outname}")
    plt.close(fig)



# =========================================================
# Main
# =========================================================
input_xyz = 'protein_nrot-5000_random.xyz'
frames = get_blocks_proteins(input_xyz)
n1, n2, n3 = 128, 363, 200

costh_list, phi_list, psi_list, cos2_list = [], [], [], []

for frame in frames:
    x = [row[1] for row in frame]
    y = [row[2] for row in frame]
    z = [row[3] for row in frame]
    # get_angles debe devolver: (costh, phi, cos2, psi)
    costh, phi, cos2, psi = get_angles(x, y, z, n1, n2, n3)
    costh_list.append(costh)
    phi_list.append(phi)
    psi_list.append(psi)
    cos2_list.append(cos2)

# Figura 3x2 (2 filas, 3 columnas)
fig, axes = plt.subplots(2, 3, figsize=(16, 8))

plot_angles_1x3(
    costh_list,   # cosθ
    phi_list,     # φ (rad)
    psi_list,     # ψ (rad)
    outname="histos_angles_1x3.png"
)


