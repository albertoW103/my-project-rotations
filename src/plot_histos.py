#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import sys


def get_blocks_proteins(input_filename):
    """
    Lee un XYZ multi-frame y devuelve una lista de frames.
    Cada frame es una lista de [label, x, y, z].
    """
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
                # Si la línea no es el número de átomos, ignora y continúa
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


def get_angles(xrot, yrot, zrot, n1, n2, n3, eps=1e-12):
    """
    Devuelve (cosθ, φ, cos²θ, ψ) para el par de vectores definidos por índices 1-based:
      - θ y φ del vector v12 = r(n1) - r(n2)
      - ψ orientado en [−π, π) usando v13 = r(n1) - r(n3) en el plano del cuerpo (conv. ZYZ)
    """
    if n1 is None or n2 is None or n3 is None:
        return None

    nu = len(xrot)
    xrot = np.asarray(xrot, float)
    yrot = np.asarray(yrot, float)
    zrot = np.asarray(zrot, float)

    # a índices 0-based
    n1 -= 1; n2 -= 1; n3 -= 1

    if not (0 <= n1 < nu and 0 <= n2 < nu and 0 <= n3 < nu):
        raise IndexError("n1/n2/n3 fuera de rango.")

    # v12
    delx = xrot[n1] - xrot[n2]
    dely = yrot[n1] - yrot[n2]
    delz = zrot[n1] - zrot[n2]

    delr2 = delx*delx + dely*dely + delz*delz
    if delr2 <= eps:
        raise ValueError("Distancia n1–n2 ~ 0: no se puede definir orientación.")

    delr = math.sqrt(delr2)
    costh = delz / delr
    # clamp por seguridad numérica
    costh = max(-1.0, min(1.0, costh))

    phi   = math.atan2(dely, delx)         # ∈ (−π, π]
    costh2 = costh*costh

    # v13 para ψ
    delx13 = xrot[n1] - xrot[n3]
    dely13 = yrot[n1] - yrot[n3]
    delz13 = zrot[n1] - zrot[n3]

    # seno/coseno de θ y φ
    sinth = math.sqrt(max(0.0, 1.0 - costh*costh))
    cosph = math.cos(phi)
    sinph = math.sin(phi)

    # En el polo (sinθ≈0), ψ es indefinido: devolvemos 0 de forma estable
    if sinth < eps:
        psi = 0.0
        return costh, phi, costh2, psi

    # Ejes del cuerpo tras Rz(φ)·Ry(θ):
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

    psi = math.atan2(py, px)               # ∈ (−π, π]
    return costh, phi, costh2, psi


def plot_angles_1x3(costheta_list, phi_list, psi_list, output_filename, bins=80):
    """
    Histos de densidad para cosθ ∈ [−1,1], φ/π ∈ [−1,1], ψ/π ∈ [−1,1].
    Una distribución uniforme ideal tiene densidad ≈ 0.5 en esos rangos.
    """
    costheta_list = np.asarray(costheta_list, float)
    phi_list      = np.asarray(phi_list, float) / np.pi   # [−1,1]
    psi_list      = np.asarray(psi_list, float) / np.pi   # [−1,1]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    dx = 0.1

    # cos(theta) en [-1, 1]
    ax = axes[0]
    ax.hist(costheta_list, bins=bins, range=(-1, 1), density=True, alpha=0.8, edgecolor='black')
    ax.axhline(0.5, linestyle='dashed', color='black')
    ax.set_title(r'Density of $\cos\theta$')
    ax.set_ylabel('Density')
    ax.set_xlabel(r'$\cos\theta$')
    ax.set_ylim(0, 1+dx)
    ax.set_xlim(-1-dx, 1+dx)
    ax.grid(True, alpha=0.3)

    # φ en [-π, π) → φ/π en [-1, 1)
    ax = axes[1]
    ax.hist(phi_list, bins=bins, range=(-1, 1), density=True, alpha=0.8, edgecolor='black')
    ax.axhline(0.5, linestyle='dashed', color='black')
    ax.set_title(r'Density of $\phi$')
    ax.set_ylabel('Density')
    ax.set_xlabel(r'$\phi/\pi$')
    ax.set_ylim(0, 1+dx)
    ax.set_xlim(-1-dx, 1+dx)
    ax.grid(True, alpha=0.3)

    # ψ en [-π, π) → ψ/π en [-1, 1)
    ax = axes[2]
    ax.hist(psi_list, bins=bins, range=(-1, 1), density=True, alpha=0.8, edgecolor='black')
    ax.axhline(0.5, linestyle='dashed', color='black')
    ax.set_title(r'Density of $\psi$')
    ax.set_ylabel('Density')
    ax.set_xlabel(r'$\psi/\pi$')
    ax.set_ylim(0, 1+dx)
    ax.set_xlim(-1-dx, 1+dx)
    ax.grid(True, alpha=0.3)

    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"[OUTPUT] Saved: {output_filename}")
    plt.show()


def parse_args():
    p = argparse.ArgumentParser(description="Plot histograms of Euler angles from a rotated XYZ trajectory.")
    p.add_argument("xyz_file",    type=str, help="Input XYZ trajectory (multi-frame)")
    p.add_argument("--n1",        type=int, default=128, help="1-based index for atom n1 (default: 128)")
    p.add_argument("--n2",        type=int, default=363, help="1-based index for atom n2 (default: 363)")
    p.add_argument("--n3",        type=int, default=200, help="1-based index for atom n3 (default: 200)")
    p.add_argument("--bins",      type=int, default=80, help="Number of histogram bins (default: 80)")
    p.add_argument("-o", "--out", type=str, default=None, help="Output image filename (.png)")
    return p.parse_args()


def main():
    args = parse_args()

    if not os.path.isfile(args.xyz_file):
        print(f"[ERROR] File not found: {args.xyz_file}")
        sys.exit(1)

    # nombre de salida
    if args.out is None:
        stem = os.path.splitext(os.path.basename(args.xyz_file))[0]
        out_png = f"histos_angles_1x3_{stem}.png"
    else:
        out_png = args.out

    print("[INPUT]  xyz_file =", args.xyz_file)
    print("[INPUT]  n1, n2, n3 =", args.n1, args.n2, args.n3)
    print("[INPUT]  bins     =", args.bins)
    print("[OUTPUT] figure   =", out_png)

    frames = get_blocks_proteins(args.xyz_file)
    if not frames:
        print("[ERROR] No frames read from XYZ.")
        sys.exit(1)

    costh_list, phi_list, psi_list, cos2_list = [], [], [], []

    for frame in frames:
        x = [row[1] for row in frame]
        y = [row[2] for row in frame]
        z = [row[3] for row in frame]
        vals = get_angles(x, y, z, args.n1, args.n2, args.n3)
        if vals is None:
            continue
        costh, phi, cos2, psi = vals
        costh_list.append(costh)
        phi_list.append(phi)
        psi_list.append(psi)
        cos2_list.append(cos2)

    if len(costh_list) == 0:
        print("[ERROR] No angles computed. Check indices n1/n2/n3.")
        sys.exit(1)

    plot_angles_1x3(costh_list, phi_list, psi_list, out_png, bins=args.bins)


if __name__ == "__main__":
    main()

