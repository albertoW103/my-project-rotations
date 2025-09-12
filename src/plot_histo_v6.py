#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import math

# ---------------------------
# Lectura de XYZ por frames
# ---------------------------
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

# ---------------------------
# Ángulos desde dos puntos
# ---------------------------
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
    if psi < 0.0:
        psi += 2.0*math.pi
    
    # definition of the fectors
    
    return costh, phi, costh2, psi

# ---------------------------
# Lectura de data.dat
# ---------------------------
#def read_data_dat(path="data.dat"):
#    thetas, phis, psis = [], [], []
#    with open(path, "r") as f:
#        for line in f:
#            if not line.strip():
#                continue
#            vals = line.split()
#            th  = float(vals[0])
#            ph  = float(vals[1])
#            ps  = float(vals[2])
#            thetas.append(th); phis.append(ph); psis.append(ps)
#    thetas = np.asarray(thetas, float)
#    phis   = np.asarray(phis, float)
#    psis   = np.asarray(psis, float)
#    # Derivados:
#    costh  = np.cos(thetas)
#    cos2th = costh**2
#    return costh, phis, cos2th

# ---------------------------
# Helper: construir edges a partir de valores discretos
# ---------------------------
#def make_edges_from_values(values, low, high):
#    """
#    A partir de una lista/array de valores 'discretos' dentro [low, high],
#    genera bordes de bins para que cada valor caiga en su propio bin.
#    """
#    v = np.unique(np.asarray(values, float))
#    v = v[(v >= low) & (v <= high)]
#    if v.size < 2:
#        # degenerado: usar un solo bin
#        return np.array([low, high], float)
#    mids = 0.5 * (v[1:] + v[:-1])  # puntos medios entre adyacentes
#    # extrapolar extremos
#    left  = v[0]  - (mids[0]   - v[0])
#    right = v[-1] + (v[-1]     - mids[-1])
#    edges = np.concatenate([[left], mids, [right]])
#    # recortar por seguridad
#    edges[0]  = max(edges[0], low)
#    edges[-1] = min(edges[-1], high)
#    return edges

# ---------------------------
# Plot helper (una fila)
# ---------------------------
# ---------------------------
# Plot helper (3x2)
# ---------------------------
def plot_angles_3x2(axes_3x2, costheta, phi, psi, cos2theta,
                    bins=80, title_prefix="",
                    align_cos_edges=False, align_phi_edges=False, align_psi_edges=False):
    # --- cos(theta) en [-1, 1] ---
    ax = axes_3x2[0, 0]
    cos_vals = np.asarray(costheta, float)
    if align_cos_edges:
        edges_cos = make_edges_from_values(cos_vals, -1.0, 1.0)
        ax.hist(cos_vals, bins=edges_cos, density=True, alpha=0.75, edgecolor='black')
    else:
        ax.hist(cos_vals, bins=bins, range=(-1, 1), density=True, alpha=0.75, edgecolor='black')
    ax.axhline(0.5, linestyle='--', color='black')  # pdf uniforme para cosθ ~ U[-1,1]
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel('cosθ'); ax.set_ylabel('Density')
    ax.set_title(f'{title_prefix}Histogram of cosθ')
    ax.grid(True)

    # --- φ/π en [0, 2] ---
    ax = axes_3x2[0, 1]
    phi_vals = (np.mod(np.asarray(phi, float), 2*np.pi)) / np.pi
    if align_phi_edges:
        edges_phi = make_edges_from_values(phi_vals, 0.0, 2.0)
        ax.hist(phi_vals, bins=edges_phi, density=True, alpha=0.75, edgecolor='black')
    else:
        ax.hist(phi_vals, bins=bins, range=(0.0, 2.0), density=True, alpha=0.75, edgecolor='black')
    ax.axhline(0.5, linestyle='--', color='black')  # uniforme en [0,2]
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel('φ/π'); ax.set_title(f'{title_prefix}Histogram of φ/π')
    ax.grid(True)

    # --- ψ/π en [0, 2] ---
    ax = axes_3x2[0, 2]
    psi_vals = (np.mod(np.asarray(psi, float), 2*np.pi)) / np.pi
    if align_psi_edges:
        edges_psi = make_edges_from_values(psi_vals, 0.0, 2.0)
        ax.hist(psi_vals, bins=edges_psi, density=True, alpha=0.75, edgecolor='black')
    else:
        ax.hist(psi_vals, bins=bins, range=(0.0, 2.0), density=True, alpha=0.75, edgecolor='black')
    ax.axhline(0.5, linestyle='--', color='black')  # uniforme en [0,2]
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel('ψ/π'); ax.set_title(f'{title_prefix}Histogram of ψ/π')
    ax.grid(True)

    # --- cos²(theta) en [0, 1] con <cos²θ> ---
    ax = axes_3x2[1, 0]
    vals_c2 = np.asarray(cos2theta, float)
    ax.hist(vals_c2, bins=bins, range=(0.0, 1.0), density=True, alpha=0.75, edgecolor='black')
    mean_cos2 = float(np.mean(vals_c2))
    ax.text(0.05, 0.95, f"<cos²θ> = {mean_cos2:.6f}",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=10, fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
    ax.set_xlabel('cos²θ'); ax.set_ylabel('Density')
    ax.set_title(f'{title_prefix}Histogram of cos²θ')
    ax.grid(True)

    # Desactivar los paneles (abajo-centro) y (abajo-derecha)
    axes_3x2[1, 1].axis('off')
    axes_3x2[1, 2].axis('off')

    # --- cos²(theta) en [0, 1] con <cos²θ> ---
    ax = axes_3x2[1, 0]
    vals_c2 = np.asarray(cos2theta, float)
    ax.hist(vals_c2, bins=bins, range=(0.0, 1.0), density=True, alpha=0.75, edgecolor='black')
    mean_cos2 = float(np.mean(vals_c2))
    ax.text(0.05, 0.95, f"<cos²θ> = {mean_cos2:.6f}",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=10, fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
    ax.set_xlabel('cos²θ'); ax.set_ylabel('Density')
    ax.set_title(f'{title_prefix}Histogram of cos²θ')
    ax.grid(True)

    # ===== Panel de texto en el espacio en blanco (abajo-centro) =====
    ax_info = axes_3x2[1, 1]
    ax_info.axis('off')

    # Promedios "centrados" para que el ideal sea 0:
    # - cosθ: promedio directo → 0 ideal
    # - φ/π: centramos restando 1 (pasa de [0,2] a [-1,1]) → 0 ideal
    # - ψ/π: idem
    N = int(np.size(costheta))
    cos_vals_mean = float(np.mean(np.asarray(costheta, float)))
    phi_pi_centered_mean = float(np.mean((np.mod(np.asarray(phi, float), 2*np.pi) / np.pi) - 1.0))
    psi_pi_centered_mean = float(np.mean((np.mod(np.asarray(psi, float), 2*np.pi) / np.pi) - 1.0))

    info_text = (
        f"N rotations = {N}\n"
        f"<cosθ>       = {cos_vals_mean:.6f}\n"
        f"<φ/π - 1>    = {phi_pi_centered_mean:.6f}\n"
        f"<ψ/π - 1>    = {psi_pi_centered_mean:.6f}\n"
        f"<cos²θ>      = {mean_cos2:.6f}"
    )

    ax_info.text(
        0.02, 0.98, info_text,
        transform=ax_info.transAxes,
        ha="left", va="top",
        fontsize=12, fontweight="bold",
        bbox=dict(facecolor="white", alpha=0.95, edgecolor="black")
    )

    # Deja apagado el panel de abajo-derecha
    axes_3x2[1, 2].axis('off')


# =========================================================
# Main
# =========================================================
input_xyz = 'protein_nrot-0_grid_lhs.xyz'
frames = get_blocks_proteins(input_xyz)
n1, n2, n3 = 128, 363, 200

costh_vec, phi_vec, psi_vec, cos2_vec = [], [], [], []

for frame in frames:
    x = [row[1] for row in frame]
    y = [row[2] for row in frame]
    z = [row[3] for row in frame]
    # get_angles debe devolver: (costh, phi, cos2, psi)
    costh, phi, cos2, psi = get_angles(x, y, z, n1, n2, n3)
    costh_vec.append(costh)
    phi_vec.append(phi)
    psi_vec.append(psi)
    cos2_vec.append(cos2)

# Figura 3x2 (2 filas, 3 columnas)
fig, axes = plt.subplots(2, 3, figsize=(16, 8))

plot_angles_3x2(
    axes,
    costh_vec,                 # cosθ
    phi_vec,                   # φ (rad)
    psi_vec,                   # ψ (rad)
    cos2_vec,                  # cos²θ
    bins=80,
    title_prefix="vectors: ",
    align_cos_edges=False, align_phi_edges=False, align_psi_edges=False
)

plt.tight_layout()
outname = "histos_3x2.png"
plt.savefig(outname, dpi=300, bbox_inches='tight')
print(f"Saved: {outname}")
plt.show()

