#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import argparse


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
    Devuelve (cosθ, φ [rad], cos²θ, ψ) del vector n1→n2 (índices 1-based).
    """
    if n1 is None or n2 is None:
        return None
    nu = len(xrot)

    xrot = np.asarray(xrot, float)
    yrot = np.asarray(yrot, float)
    zrot = np.asarray(zrot, float)

    # indices a 0-based
    n1 -= 1; n2 -= 1; n3 -= 1
    
    if not (0 <= n1 < nu and 0 <= n2 < nu and 0 <= n3 < nu):
        raise IndexError("n1/n2/n3 out of range.")

    delx = xrot[n1] - xrot[n2]
    dely = yrot[n1] - yrot[n2]
    delz = zrot[n1] - zrot[n2]
    
    delr = math.sqrt(delx*delx + dely*dely + delz*dely)
    
    costh = delz / delr
    phi   = math.atan2(dely, delx)
    costh2 = costh*costh

    # Vector n1→n3
    delx13 = xrot[n1] - xrot[n3]
    dely13 = yrot[n1] - yrot[n3]
    delz13 = zrot[n1] - zrot[n3]
    
    sinth = math.sqrt(max(0.0, 1.0 - costh*costh))
    cosph = math.cos(phi)
    sinph = math.sin(phi)
    
    exx = cosph * costh
    exy = sinph * costh
    exz = -sinth
    
    eyx = -sinph
    eyy =  cosph
    eyz =  0.0
    
    px = delx13*exx + dely13*exy + delz13*exz
    py = delx13*eyx + dely13*eyy + delz13*eyz
    
    psi = math.atan2(py, px)
           
    return costh, phi, costh2, psi


def plot_angles_1x3(costheta_list, phi_list, psi_list, output_filename):
    costheta_list = np.array(costheta_list)
    phi_list      = np.array(phi_list) / np.pi
    psi_list      = np.array(psi_list) / np.pi
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    bins = 80
    dx = 0.1
    
    # cos(theta)
    ax = axes[0]
    ax.hist(costheta_list, bins=bins, range=(-1, 1), density=True, alpha=0.8, edgecolor='black')
    ax.axhline(0.5, linestyle='dashed', color='black')
    ax.set_title('Density of $\cos\\theta$')
    ax.set_ylabel('Density')
    ax.set_xlabel(r'$\cos\theta$')
    ax.set_ylim(0, 1+dx)
    ax.set_xlim(-1-dx, 1+dx)
    ax.grid(True, alpha=0.3)
    
    # phi
    ax = axes[1]
    ax.hist(phi_list, bins=bins, range=(-1, 1), density=True, alpha=0.8, edgecolor='black')
    ax.axhline(0.5, linestyle='dashed', color='black')
    ax.set_title('Density of $\phi$')
    ax.set_ylabel('Density')
    ax.set_xlabel(r'$\phi$/$\pi$')
    ax.set_ylim(0, 1+dx)
    ax.set_xlim(-1-dx, 1+dx)
    ax.grid(True, alpha=0.3)
    
    # psi
    ax = axes[2]
    ax.hist(psi_list, bins=bins, range=(-1, 1), density=True, alpha=0.8, edgecolor='black')
    ax.axhline(0.5, linestyle='dashed', color='black')
    ax.set_title('Density of $\psi$')
    ax.set_ylabel('Density')
    ax.set_xlabel(r'$\psi$/$\pi$')
    ax.set_ylim(0, 1+dx)
    ax.set_xlim(-1-dx, 1+dx)
    ax.grid(True, alpha=0.3)

    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_filename}")
    plt.show()


# =========================================================
# Main
# =========================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot histograms of Euler angles from rotated XYZ.")
    parser.add_argument("input_file", type=str, help="Input multi-frame XYZ file (e.g., 4F5S_nrot-1000_random.xyz)")
    parser.add_argument("-n1", type=int, default=128, help="Index of atom 1 (1-based)")
    parser.add_argument("-n2", type=int, default=363, help="Index of atom 2 (1-based)")
    parser.add_argument("-n3", type=int, default=200, help="Index of atom 3 (1-based)")
    args = parser.parse_args()

    frames = get_blocks_proteins(args.input_file)

    costh_list, phi_list, psi_list, cos2_list = [], [], [], []
    for frame in frames:
        x = [row[1] for row in frame]
        y = [row[2] for row in frame]
        z = [row[3] for row in frame]
        costh, phi, cos2, psi = get_angles(x, y, z, args.n1, args.n2, args.n3)
        costh_list.append(costh)
        phi_list.append(phi)
        psi_list.append(psi)
        cos2_list.append(cos2)

    output_filename = "histos_angles_1x3.png"
    plot_angles_1x3(costh_list, phi_list, psi_list, output_filename)

