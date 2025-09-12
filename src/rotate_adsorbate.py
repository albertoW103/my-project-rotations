#!/usr/bin/env python3


"""
rotate_adsorbate.py

Generate rotated configurations of a molecule/protein using Euler angles (ZYZ). 
Supports multiple sampling modes:
- grid      : open intervals, midpoints
- grid_2    : endpoints included
- grid_lhs  : Latin Hypercube stratified sampling (exact nrot)
- random    : Haar-like random sampling

Outputs:
- An XYZ trajectory file with rotated coordinates.
- A 'data.dat' file with the angles (theta, phi, psi) used.

Dependencies:
- numpy, scipy, matplotlib

Example usage:
    python rotate_adsorbate.py --input protein.xyz --mode grid_lhs --nrot 100
"""



import math
import numpy as np
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
import argparse
import sys
import os

# ---------------------------
# Utilities
# ---------------------------

def ask_positive_int(prompt: str) -> int:
    """
    Keep asking until the user provides a positive integer (>0).
    """
    while True:
        try:
            value_str = input(prompt).strip()
            value = int(value_str)
            if value <= 0:
                raise ValueError("Value must be a positive integer greater than zero.")
            return value
        except ValueError as e:
            print(f"[Error]: {e}. Please try again.")

def ask_mode() -> str:
    """
    Keep asking until the user provides a valid mode.
    Valid options: 'grid' or 'random'.
    """
    while True:
        mode = input("Mode of rotation [grid / random]: ").strip().lower()
        if mode in ["random", "grid"]:
            return mode
        else:
            print("[Error]: You need to specify a valid mode (grid or random). Please try again.")

def count_unique_angles(angles, decimals=None):
    """
    Cuenta valores únicos en `angles`.
    - decimals=None: unicidad exacta de float (prácticamente será nrot).
    - decimals=k: redondea a k decimales antes de contar (útil si querés tolerancia).
    """
    arr = np.asarray(angles, dtype=float)
    if decimals is not None:
        arr = np.round(arr, decimals)
    return np.unique(arr).size


def plot_three_angles(theta_list, phi_list, psi_list, radius=2*np.pi, filename=None):
    """
    Una sola figura con 3 subplots (θ, φ, ψ) dibujados como flechas en el círculo.
    """
    def draw(ax, angles, title):
        ang = np.asarray(angles, float) % (2*np.pi)
        t = np.linspace(0, 2*np.pi, 400)
        xc, yc = radius*np.cos(t), radius*np.sin(t)
        x = radius*np.cos(ang); y = radius*np.sin(ang)
        ax.plot(xc, yc, lw=1.0, color='black')
        ax.axhline(0, lw=0.5, color='gray'); ax.axvline(0, lw=0.5, color='gray')
        for xi, yi in zip(x, y):
            ax.arrow(0, 0, xi, yi, length_includes_head=True,
                     head_width=0.12, head_length=0.28, linewidth=0.9)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(-1.1*radius, 1.1*radius); ax.set_ylim(-1.1*radius, 1.1*radius)
        ax.set_title(title); ax.set_xticks([]); ax.set_yticks([]); ax.grid(True, alpha=0.3, lw=0.4)

    fig, axs = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    draw(axs[0], theta_list, "θ (theta)")
    draw(axs[1], phi_list,   "φ (phi)")
    draw(axs[2], psi_list,   "ψ (psi)")

    if filename:
        fig.savefig(filename, dpi=200, bbox_inches="tight")
    plt.close(fig)
    

def get_midpoints(rmin, rmax, dimr):
    """
    Devuelve (delta, midpoints) para 'dimr' bins iguales en el intervalo ABIERTO (rmin, rmax),
    o sea, centros que evitan los extremos.

    Comportamiento:
      - dimr <= 0 : (None, array([]))
      - dimr == 1 : (None, array([ (rmin+rmax)/2 ]))  -> sin rango, un solo punto medio
      - dimr >= 2 : (delta, array de midpoints)
    """
    if dimr is None or dimr <= 0:
        return None, np.array([], dtype=float)
    if dimr == 1:
        mid = 0.5 * (rmin + rmax)
        return None, np.array([mid], dtype=float)
    delta = (rmax - rmin) / dimr
    midpoints = rmin + (np.arange(dimr) + 0.5) * delta
    return delta, midpoints


def print_delta(label, delta, arr, denom=1.0, suffix=""):
    """
    Imprime delta, o el punto medio si no hay rango (delta=None).
    """
    if delta is None:
        if arr.size:
            val = float(arr[0]) / denom
            print(f"{label} = no range (midpoint = {val:.2f}{suffix})")
        else:
            print(f"{label} = no range (empty)")
    else:
        print(f"{label} = {delta/denom:.2f}{suffix}")
        

def rota(xsv: np.ndarray,
	 ysv: np.ndarray,
	 zsv: np.ndarray,
	 theta, phi, psi):
    """
    Rotación con Euler ZYZ. Esta función RECIBE (theta, phi, psi) y aplica:
      R = Rz(phi) · Ry(theta) · Rz(psi)
    nota: search convention ZY′Z″:
    se rota primero alrededor de Z, luego alrededor de Y′ (ya rotado), y por último alrededor de Z″ (el Z después de las dos primeras).
    ver este video: https://www.youtube.com/watch?v=N7AVc5yYX-k
    el primer angulo seria phi, el segundo theta, y el tercero psi.
    
    Parameters:
        x (np.ndarray): X coordinates, shape (N,).
        y (np.ndarray): Y coordinates, shape (N,).
        z (np.ndarray): Z coordinates, shape (N,).
        theta (float): Polar angle θ in radians.
        phi (float): Azimuthal angle φ in radians.
        psi (float): Roll angle ψ in radians.
    Returns:
        xrot (np.ndarray): Rotated X, shape (N,).
        yrot (np.ndarray): Rotated Y, shape (N,).
        zrot (np.ndarray): Rotated Z, shape (N,).
    
    """
    coords = np.column_stack((xsv, ysv, zsv))                              # (N,3)
    rotation = R.from_euler("ZYZ", [phi, theta, psi], degrees=False)       # φ, θ, ψ en rad
    coords_rot = rotation.apply(coords)                                    # get xyz rotated coordenates
    xrot, yrot, zrot = coords_rot[:,0], coords_rot[:,1], coords_rot[:,2]   # split coordenates
    return  xrot, yrot, zrot                                               # np.arrays


def adsorbate_rot(
        xyz_output_filename,
        nu, header, res,
        x:np.ndarray,
        y:np.ndarray,
        z:np.ndarray,
        ntheta=None, nphi=None, npsi=None,
        mode=None,
        nrot=None
    ):
    """
    Generate rotations and write XYZ frames. Also save (theta, phi, psi) to 'data.dat'.

    Parameters
    ----------
    xyz_output_filename : str
        Output XYZ file path where all rotated frames will be written.
    nu : int
        Number of points/atoms.
    header : str
        Header line to write in each XYZ frame.
    res : list[str]
        Labels per point (length nu).
    x, y, z : np.ndarray
        Coordinate arrays of shape (nu,).
    ntheta, nphi, npsi : int or None
        Grid sizes for theta, phi, psi (used in 'grid' / 'grid_2' modes).
    mode : {"grid", "grid_2", "random"}
        Rotation generation mode.
    nrot : int or None
        Number of random rotations (only used when mode == "random").
    """
    # defs:
    pi = math.pi    
    rng = np.random.default_rng(seed=0) # random number generator
    
    ####################################
    # Centered coordinates (one-time shift)
    ####################################
    # center of mass:
    xcm = ycm = zcm = 0.0
    for iu in range(nu):
        xcm += x[iu]
        ycm += y[iu]
        zcm += z[iu]
    xcm /= nu
    ycm /= nu
    zcm /= nu
    
    # shift to COM (labeled as shifted vector, sv):
    xsv = x - xcm   # arrays
    ysv = y - ycm   # arrays
    zsv = z - zcm   # arrays
    
    #########################################
    #  Get angles to rotate the molecule
    #########################################
    
    if mode == "grid":
        ###################################################
        # Sampling with midpoints (open intervals):
        # cos(theta) ∈ (-1, 1),  phi ∈ (0, 2π),  psi ∈ (0, 2π)
        ###################################################
        print('method used: grid')
        
        # require all grid sizes:
        if None in (ntheta, nphi, npsi):
            print("grid: missing parameters (ntheta/nphi/npsi). Doing nothing.")
            return
        
        # warning if nx are < 1:
        if ntheta <= 0 or nphi <= 0 or npsi <= 0:
            print("grid: all ntheta/nphi/npsi must be > 0.")
            return
            
        # warning if theta sampling is too coarse:
        if ntheta < 10:
            print("grid: warning — coarse sampling in cos(theta); results may be biased.")
        
        # empty list to store angles:
        triplet_list = []
    
        # midpoints (open intervals):
        delta_costheta, costheta_arr = get_midpoints(-1.0, 1.0, ntheta)  # cosθ -> (-1,1)
        delta_phi, phi_arr           = get_midpoints( 0.0, 2*pi, nphi)   # φ -> (0,2π)
        delta_psi, psi_arr           = get_midpoints( 0.0, 2*pi, npsi)   # ψ -> (0,2π)
        
        # examples:
        # Δθ = 0.5,  midpoints: [-0.75 -0.25  0.25  0.75]
        # Δφ =  π/2, midpoints: [π/4, 3π/4, 5π/4, 7π/4]
        # Δψ =  π/2, midpoints: [π/4, 3π/4, 5π/4, 7π/4]
        
        #psi_arr = psi_arr + delta_psi/2
        
        # Note: this is Δ(cosθ), not Δθ:
        print_delta("Δcosθ", delta_costheta, costheta_arr)
        print_delta("Δφ",    delta_phi,      phi_arr, math.pi, " π")
        print_delta("Δψ",    delta_psi,      psi_arr, math.pi, " π")
        
        # to lists:
        theta_list = np.arccos(costheta_arr).tolist()
        phi_list = phi_arr.tolist()
        psi_list = psi_arr.tolist()
        
        # build triplets:
        for theta in theta_list:
            for phi in phi_list:
                for psi in psi_list:
                    triplet_list.append((theta, phi, psi))
        
        length = len(triplet_list)
        print(length, 'rotations will be generated')

        unique_theta_exact = count_unique_angles(theta_list, decimals=None)
        unique_phi_exact   = count_unique_angles(phi_list, decimals=None)
        unique_psi_exact   = count_unique_angles(psi_list, decimals=None)
        print(f"unique θ count: {unique_theta_exact}")
        print(f"unique φ count: {unique_phi_exact}")
        print(f"unique ψ count: {unique_psi_exact}")
                
        # visualize angles:
        plot_three_angles(theta_list, phi_list, psi_list, filename=f"angles_three_{mode}.png")

    elif mode == "grid_2":
        ###################################################
        # Endpoints included:
        # cosθ ∈ [-1, 1] (includes ±1)  → θ = arccos(cosθ)
        # φ, ψ ∈ [0, 2π) (include 0, exclude 2π to avoid duplicates)
        ###################################################
        print('method used: grid_2 (with endpoints)')

        if None in (ntheta, nphi, npsi):
            print("grid_2: missing parameters (ntheta/nphi/npsi). Doing nothing.")
            return
        
        if ntheta <= 0 or nphi <= 0 or npsi <= 0:
            print("grid_2: all ntheta/nphi/npsi must be > 0.")
            return

        # warning if theta sampling is too coarse:
        if ntheta < 10:
            print("grid_2: warning — coarse sampling in cos(theta); results may be biased.")
            
        # cosθ with endpoints; φ and ψ with 0 included and 2π excluded:
        if ntheta == 1:
            costheta_arr = np.array([0.0])                  # un solo valor: θ = π/2
            delta_costheta = None
        else:
            costheta_arr = np.linspace(-1.0, 1.0, ntheta, endpoint=True)
            delta_costheta = 2.0 / (ntheta - 1)

        phi_arr = np.linspace(0.0, 2*pi, nphi, endpoint=False)
        psi_arr = np.linspace(0.0, 2*pi, npsi, endpoint=False)
        delta_phi = (2*pi) / nphi if nphi > 0 else None
        delta_psi = (2*pi) / npsi if npsi > 0 else None

        # grid info:
        print_delta("Δcosθ", delta_costheta, costheta_arr)
        print_delta("Δφ",    delta_phi,      phi_arr, math.pi, " π")
        print_delta("Δψ",    delta_psi,      psi_arr, math.pi, " π")

        # to lists:
        theta_list = np.arccos(np.clip(costheta_arr, -1.0, 1.0)).tolist()
        phi_list   = phi_arr.tolist()
        psi_list   = psi_arr.tolist()

        # build triplets:
        triplet_list = []
        for theta in theta_list:
            for phi in phi_list:
                for psi in psi_list:
                    triplet_list.append((theta, phi, psi))  # triple of angles to generate rotations
        
        length = len(triplet_list)
        print(length, 'rotations will be generated (grid_2)')

        # visualize angles:
        plot_three_angles(theta_list, phi_list, psi_list, filename=f"angles_three_{mode}.png")
                
    elif mode == 'random':
        ##########################################################
        # Haar-like sampling via independent ZYZ parameters:
        # θ ∈ [0, π),  φ ∈ [0, 2π),  ψ ∈ [0, 2π)
        ##########################################################
        print('method used: random')
        
        # warning if nrot < 1:
        if (nrot is None) or (nrot < 1):
            print("random: 'nrot' must be > 0. Doing nothing.")
            return
        
        # empty list to store angles:
        triplet_list = []
        
        for _ in range(nrot):
            # get random numbers ro defined angles:
            u1, u2, u3 = rng.random(), rng.random(), rng.random()   # rng.random() in [0,1)
            
            # theta:
            theta = np.arccos(2.0*u1 - 1.0)   # get cosθ ~ U[-1,1)
            
            # phi and psi (the same for both):
            phi   = 2.0*pi*u2                 # U[0,2π)
            psi   = 2.0*pi*u3                 # U[0,2π)
            
            # append:
            triplet_list.append((theta, phi, psi))
            
        length = len(triplet_list)
        print(length, 'rotations will be generated')
        
        # plot:
        theta_list = [t[0] for t in triplet_list]
        phi_list   = [t[1] for t in triplet_list]
        psi_list   = [t[2] for t in triplet_list]
        
        unique_theta_exact = count_unique_angles(theta_list, decimals=None)
        unique_phi_exact   = count_unique_angles(phi_list, decimals=None)
        unique_psi_exact   = count_unique_angles(psi_list, decimals=None)
        print(f"unique θ count: {unique_theta_exact}")
        print(f"unique φ count: {unique_phi_exact}")
        print(f"unique ψ count: {unique_psi_exact}")
                
        plot_three_angles(theta_list, phi_list, psi_list, filename=f"angles_three_{mode}.png")
        
    else:
        raise ValueError("mode must be one of {'grid', 'grid_2', 'random'}.")
    
    #########################################################
    # Apply rotations (one per triplet) and write XYZ
    #########################################################
    # ope file:
    with open(xyz_output_filename, 'w') as file:
        irot = 0
        for theta, phi, psi in triplet_list:   # the rotation is made from triplet
            xrot, yrot, zrot = rota(xsv, ysv, zsv, theta, phi, psi)
            
            # write molecule
            file.write(f"{nu}\n")
            file.write(f"{header}\n")
            for iu in range(nu):
                file.write(f"{res[iu]:<8s}"
                           f"{xrot[iu]:10.4f}"
                           f"{yrot[iu]:10.4f}"
                           f"{zrot[iu]:10.4f}\n")
            irot += 1
    print('n rotations =', irot, ' are generated')
    
    #########################################################
    # Save angles used for plotting on the sphere
    #########################################################
    outfile = "data.dat"
    with open(outfile, "w") as fa:
        for theta, phi, psi in triplet_list:  # the rotation is made from triplet
            fa.write(f"{theta:.10f}  {phi:.10f}  {psi:.10f}\n")
    print("Triplet angles saved as data.dat")
    

def read_xyz_file(input_filename):
    """
    Read a simple XYZ (N, header, then label x y z) and return (nu, header, res, x, y, z).
    
    Parameters:
        input_filename (str): Path to input .xyz file.
    Returns:
        nu (int): Number of points/atoms.
        header (str): XYZ header line.
        res (list[str]): Labels per point, length nu.
        x (np.ndarray): X coordinates, shape (nu,).
        y (np.ndarray): Y coordinates, shape (nu,).
        z (np.ndarray): Z coordinates, shape (nu,).
        
    """
    
    # Open file:
    with open(input_filename, 'r') as file:
        # Number of atoms/points:
        nu = int(file.readline())
        
        # Header line:
        header = file.readline().rstrip("\n")
        
        # Accumulators:
        res, x, y, z = [], [], [], []
        
        # Read N coordinate lines:
        for _ in range(nu):
            line = file.readline().split()
            resi = line[0]
            xi = float(line[1])
            yi = float(line[2])
            zi = float(line[3])
            
            res.append(resi)
            x.append(xi)
            y.append(yi)
            z.append(zi)
    
    # Convert to numpy arrays:
    x = np.array(x)   # arrays
    y = np.array(y)   # arrays
    z = np.array(z)   # arrays
    
    return nu, header, res, x, y, z


def angles_to_sphere_points():
    '''
    Plot a unit sphere and overlay points computed from the spherical angles
    theta (polar) and phi (azimuthal) read from 'data.dat'.
    '''
    print('Plotting points using theta and phi angles')
    
    # Extract angles from file:
    input_filename="data.dat"
    thetas, phis = [], []
    with open(input_filename, "r") as file:
        for line in file:
            vals = line.split()           
            theta = float(vals[0])
            phi = float(vals[1])
            thetas.append(theta)
            phis.append(phi)
    
    # Convert to arrays:
    th = np.asarray(thetas)
    ph = np.asarray(phis)
    
    # Convert spherical (theta, phi) to Cartesian on the unit sphere:
    x = np.sin(th)*np.cos(ph)
    y = np.sin(th)*np.sin(ph)
    z = np.cos(th)
    
    # Initialize the plot:
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")
    
    # Reference sphere mesh:
    u = np.linspace(0, 2*np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, linewidth=0.3, alpha=0.2)
    
    # Points:
    ax.scatter(x, y, z, s=12, alpha=0.85)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlim(-1,1); ax.set_ylim(-1,1); ax.set_zlim(-1,1)
    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    ax.set_title("Points on sphere (θ, φ)")
    plt.show()


########################################################
########################################################

input_filename = input("Nombre del archivo XYZ (ej: protein.xyz): ").strip()

if not os.path.isfile(input_filename):
    print(f"[Error]: file '{input_filename}' was not found.")
    sys.exit(1)
    
nu, header, res, x, y, z = read_xyz_file(input_filename)

mode = ask_mode()
if mode == "random":
    nrot_str             = input("Numer of rotations: ").strip()
    nrot                 = int(nrot_str)
    ntheta = nphi = npsi = None
    xyz_output_filename  = f"{input_filename.split('.')[0]}_nrot-{nrot}_{mode}.xyz"
    
elif mode == "grid":
    ntheta = ask_positive_int("Number of divisions in θ: ")
    nphi   = ask_positive_int("Number of divisions in φ: ")
    npsi   = ask_positive_int("Number of divisions in ψ: ")
    
    nrot = ntheta * nphi * npsi
    xyz_output_filename = f"{input_filename.split('.')[0]}_nrot-{nrot}_{mode}.xyz"
    
# get rotations:
print(f"name of the protein: {input_filename.split('.')[0]}")
nrotx = adsorbate_rot(xyz_output_filename, nu, header, res, x, y, z, ntheta, nphi, npsi, mode, nrot)

# plot:
angles_to_sphere_points()


########################################################
# un esquema tipico en grilla es:
# ntheta
# nphi = 2*ntheta
# npsi = 2*ntheta
# 6 x 12 x 12 = 864 parecen ser las rotaciones de stefy
########################################################





exit()





