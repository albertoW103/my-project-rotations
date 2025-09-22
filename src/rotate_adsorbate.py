#!/usr/bin/env python3


"""
rotate_adsorbate.py

Generate rotated configurations of a molecule/protein using ZYZ Euler angles.
Sampling modes:
- grid      : open intervals (midpoints)
- grid_2    : endpoints included
- random    : Haar-like random sampling

Outputs:
- XYZ trajectory file with rotated coordinates
- 'data.dat' with the angles used (theta, phi, psi)

Dependencies:
- numpy, scipy, matplotlib

Example:
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


def get_cli_args():
    parser = argparse.ArgumentParser(description="Generate rotated configurations of a molecule/protein")
    parser.add_argument("input_file", type=str, help="XYZ input file (e.g. protein.xyz)")
    parser.add_argument("-mode", type=str, choices=["grid", "random", "kuffner"], required=True,
                        help="Rotation mode: grid, random, or kuffner")
    parser.add_argument("-nrot", type=int, default=None,
                        help="Number of rotations (for random/kuffner)")
    parser.add_argument("-ntheta", type=int, default=None, help="Divisions in θ (grid only)")
    parser.add_argument("-phi",    type=int, default=None, help="Divisions in φ (grid only)")
    parser.add_argument("-psi",    type=int, default=None, help="Divisions in ψ (grid only)")
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        print(f"[Error]: file '{args.input_file}' was not found.")
        sys.exit(1)

    if args.mode in ("random", "kuffner"):
        if args.nrot is None or args.nrot < 1:
            print("[Error]: you must provide a positive -nrot for random/kuffner modes.")
            sys.exit(1)

    if args.mode == "grid":
        if args.ntheta is None or args.phi is None or args.psi is None:
            print("[Error]: you must provide -ntheta, -phi and -psi for grid mode.")
            sys.exit(1)

    return args


def check_unique_rotations(triplet_list, convention=None, tol=1e-12):
    """
    Detecta rotaciones repetidas comparando matrices de rotación cuantizadas.
    'triplet_list' SIEMPRE viene como (theta, phi, psi) en radianes.
    - Si convention == "ZYZ", SciPy espera (phi, theta, psi)  -> reordenar
    - Si convention == "ZYX", SciPy espera (yaw, pitch, roll) -> usar tal cual (theta, phi, psi)
    """
    if not triplet_list:
        return 0, 0, 0, np.array([], dtype=int)

    tps = np.asarray(triplet_list, dtype=float)  # (N, 3) con columnas (theta, phi, psi)

    if convention == "ZYZ":
        # SciPy espera (phi, theta, psi)
        angles = np.column_stack([tps[:,1], tps[:,0], tps[:,2]])
    elif convention == "ZYX":
        # SciPy espera (yaw, pitch, roll) = (theta, phi, psi) aquí
        angles = tps
    else:
        raise ValueError("Unsupported Euler sequence (expected 'ZYZ' or 'ZYX').")

    Rm = R.from_euler(convention, angles, degrees=False).as_matrix().reshape(len(tps), 9)
    flat_quant = np.round(Rm / tol).astype(np.int64)
    _, keep_idx, _ = np.unique(flat_quant, axis=0, return_index=True, return_inverse=True)

    n_total  = len(tps)
    n_unique = len(keep_idx)
    n_dupes  = n_total - n_unique
    return n_total, n_unique, n_dupes, np.sort(keep_idx)


def ask_positive_int(prompt: str) -> int:
    """
    Prompt until the user provides a positive integer (> 0).
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
    while True:
        mode = input("Mode of rotation [grid / random / kuffner]: ").strip().lower()
        if mode in ["random", "grid", "kuffner"]:
            return mode
        else:
            print("[Error]: valid modes are: grid, random, kuffner. Please try again.")


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


def plot_three_angles(theta_list, phi_list, psi_list, filename):
    """
    Create a single figure with 3 subplots (θ, φ, ψ).
    Each subplot shows the input angles as arrows drawn from the origin to the circumference of a circle.
    """

    radius = 2 * np.pi  # radius of the reference circle (arbitrary choice, only affects arrow length)

    def draw(ax, angles, title):
        # Convert input list to numpy array and normalize angles to [0, 2π)
        ang = np.asarray(angles, float) % (2 * np.pi)

        # Generate the circle outline for reference
        t = np.linspace(0, 2 * np.pi, 400)
        xc, yc = radius * np.cos(t), radius * np.sin(t)

        # Convert each angle to its (x, y) endpoint on the circle
        x = radius * np.cos(ang)
        y = radius * np.sin(ang)

        # Plot the circle outline
        ax.plot(xc, yc, lw=1.0, color='black')

        # Add horizontal and vertical axes for orientation
        ax.axhline(0, lw=0.5, color='gray')
        ax.axvline(0, lw=0.5, color='gray')

        # Draw one arrow for each angle, from the origin to its (x, y) endpoint
        for xi, yi in zip(x, y):
            ax.arrow(
                0, 0, xi, yi,
                length_includes_head=True,  # include arrow head in the total length
                head_width=0.12,            # width of the arrow head
                head_length=0.28,           # length of the arrow head
                linewidth=0.9               # line thickness
            )

        # Keep aspect ratio square (circle not distorted)
        ax.set_aspect('equal', adjustable='box')

        # Set limits slightly larger than radius so arrows fit nicely
        ax.set_xlim(-1.1 * radius, 1.1 * radius)
        ax.set_ylim(-1.1 * radius, 1.1 * radius)

        # Add title for this subplot and clean up axes
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(True, alpha=0.3, lw=0.4)  # light background grid

    # Create figure with 3 subplots (side by side)
    fig, axs = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    # Draw each set of angles on its corresponding subplot
    draw(axs[0], theta_list, "θ (theta)")
    draw(axs[1], phi_list,   "φ (phi)")
    draw(axs[2], psi_list,   "ψ (psi)")

    # Save the figure to file and close it
    fig.savefig(filename, dpi=300)
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


def angles_to_sphere_points_general(input_filename, mode, output_filename):
    """
    Dibuja en la esfera la orientación del eje Z del cuerpo (R * ez),
    usando la convención correcta según el modo:
      - 'random' y 'grid'  -> ZYZ  con orden (phi, theta, psi)
      - 'kuffner'          -> ZYX  con orden (theta, phi, psi)  (yaw, pitch, roll)
    """
    # 1) leer ángulos
    # Extract angles from file:
    thetas, phis, psis = [], [], []
    with open(input_filename, "r") as file:
        for line in file:
            # split:
            vals = line.split()       
            
            # get values:
            theta = float(vals[0])
            phi = float(vals[1])
            psi = float(vals[2])
            
            # append:
            thetas.append(theta)
            phis.append(phi)
            psis.append(psi)
    
    # Convert to arrays:
    thetas = np.array(thetas)
    phis = np.array(phis)
    psis = np.array(psis)
    
    # 2) construir rotaciones segun convención:
    if mode in ("random","grid"):
        convention = "ZYZ"
        # SciPy's ZYZ expects (φ, θ, ψ)
        # array de N tripletes
        euler = np.column_stack([phis, thetas, psis])  # ZYZ espera (phi, theta, psi)
        
    elif mode == "kuffner":
        convention = "ZYX"
        # SciPy's ZYX expects (yaw, pitch, roll) = (θ, φ, ψ) here
        # array de N tripletes
        euler = np.column_stack([thetas, phis, psis])  # ZYX espera (yaw, pitch, roll)
        
    else:
        raise ValueError("mode must be 'random', 'grid', or 'kuffner'")
    
    # construye las rotaciones en 3D a partir de tus ángulos de Euler usando SciPy.
    Rm = R.from_euler(convention, euler, degrees=False)
    
    # ---- 3) Apply rotations to e_z to get the oriented body Z-axis ----
    ez = np.array([[0.0, 0.0, 1.0]])  # shape (1,3); apply() broadcasts to N×3
    pts = Rm.apply(ez)                # shape (N,3)
    x, y, z = pts[:,0], pts[:,1], pts[:,2]

    # 4) graficar
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection="3d")

    u = np.linspace(0, 2*np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, linewidth=0.3, alpha=0.2)

    ax.scatter(x, y, z, s=12, alpha=0.85)
    ax.set_box_aspect((1,1,1))
    ax.set_xlim(-1,1); ax.set_ylim(-1,1); ax.set_zlim(-1,1)
    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    ax.set_title(f"Points on sphere (mode = {mode}, seq = {convention})")
    fig.savefig(output_filename, dpi=300, bbox_inches="tight")
    plt.show()
    
def rota(xsv: np.ndarray,
	 ysv: np.ndarray,
	 zsv: np.ndarray,
	 theta, phi, psi,
	 mode):  # 'the most used: ZYZ'
    """
    Aplica una rotación a las coordenadas (xsv, ysv, zsv) usando ángulos de Euler.
    - Para mode in {'random','grid'} usa convención ZYZ (orden intrínseco): Rz(phi) · Ry(theta) · Rz(psi)
    - Para mode == 'kuffner' usa convención ZYX (RPY/yaw–pitch–roll): Rz(theta) · Ry(phi) · Rx(psi)

    Devuelve arrays rotados xrot, yrot, zrot.
    """
    coords = np.column_stack((xsv, ysv, zsv))                              # (N,3)
    if mode in ('random', 'grid'):
        convention = 'ZYZ'
        # SciPy espera el orden de ángulos siguiendo la secuencia:
        euler = [phi, theta, psi] # ZYZ → (phi, theta, psi)
    elif mode == 'kuffner':
        convention = 'ZYX'
        euler = [theta, phi, psi]  # ZYX → (yaw=theta, pitch=phi, roll=psi)
    else:
        raise ValueError("Unsupported mode (expected 'random', 'grid', or 'kuffner').")
        
    rotation = R.from_euler(convention, euler, degrees=False)  # φ, θ, ψ en rad
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
        Grid sizes for theta, phi, psi (used in 'grid').
    mode : {"grid",  "random"}
        Rotation generation mode.
    nrot : int or None
        Number of random rotations (only used when mode == "random").
    """
    # defs:
    pi = math.pi    
        
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
        print('\n............................................')
        print('method used: grid (midpoints, open intervals)')
        print(f"grid: divisions = ntheta={ntheta}, nphi={nphi}, npsi={npsi}")
        print('............................................\n')
        
        # here, we defined the convention to use:
        convention = "ZYZ"
        
        # require all grid sizes:
        if None in (ntheta, nphi, npsi):
            print("[grid][abort]: missing nθ/nφ/nψ")
            return
        
        # warning if nx are < 1:
        if ntheta <= 0 or nphi <= 0 or npsi <= 0:
            print("[grid][abort]: all nθ/nφ/nψ must be > 0")
            return
            
        # warning if theta sampling is too coarse:
        if ntheta < 10:
            print("[grid][warn ]: coarse sampling in cosθ; results may be biased")
        
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
        
    elif mode == 'random':
        ##########################################################
        # Haar-like sampling via independent ZYZ parameters:
        # θ ∈ [0, π),  φ ∈ [0, 2π),  ψ ∈ [0, 2π)
        ##########################################################
        print('\n............................................')
        print('method used: random (Haar-like sampling)')
        print('............................................\n')

        # here, we defined the convention to use:
        convention = "ZYZ"
                
        seed = 0
        rng  = np.random.default_rng(seed) # random number generator
        print(f"random: RNG seed = {seed}")
        
        # warning if nrot < 1:
        if (nrot is None) or (nrot < 1):
            print("[random][abort]: 'nrot' must be > 0")
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
            
    elif mode == 'kuffner':
        ##########################################################
        # Kuffner (2004) algorithm for uniform random Euler angles
        # Roll–Pitch–Yaw convention (θ, φ, η)
        ##########################################################
        print('\n............................................')
        print('method used: kuffner (uniform RPY random sampling)')
        print('............................................\n')

        # here, we defined the convention to use:
        convention = "ZYX"
                
        seed = 0
        rng  = np.random.default_rng(seed)  # random number generator
        print(f"kuffner: RNG seed = {seed}")
        
        if (nrot is None) or (nrot < 1):
            print("[kuffner][abort]: 'nrot' must be > 0")
            return
        
        triplet_list = []
        
        for _ in range(nrot):
            # Generate random numbers:
            u1, u2, u3, u4 = rng.random(), rng.random(), rng.random(), rng.random()
            
            # Algorithm from pseudocode:
            theta  = 2.0*pi*u1 - pi                        # [-pi, pi)
            phi    = np.arccos(1.0 - 2.0*u2) + (pi/2.0)    # [-pi/2, pi/2)
            
            # reflection adjustment:
            if u3 < 0.5:
                if phi < pi:
                    phi = phi + pi
                else:
                    phi = phi - pi
            
            psi = 2.0*pi*u4 - pi                          # [-pi, pi)
            
            # store triplet:
            triplet_list.append((theta, phi, psi))
            
    else:
        raise ValueError("mode must be one of {'grid', 'grid_2', 'random'}.")
    
    
    ###################################
    #
    ###################################
    theta_list = [t[0] for t in triplet_list]
    phi_list   = [t[1] for t in triplet_list]
    psi_list   = [t[2] for t in triplet_list]
    
        
    n_total, n_unique, n_dupes, keep_idx = check_unique_rotations(triplet_list, convention)
    print(f"[dedup] total={n_total}, unique={n_unique}, duplicated={n_dupes}")
    
    # calculate lenght:
    length = len(triplet_list)
    unique_theta_exact = count_unique_angles(theta_list, decimals=None)
    unique_phi_exact   = count_unique_angles(phi_list, decimals=None)
    unique_psi_exact   = count_unique_angles(psi_list, decimals=None)
    
    print(f"total rotations to generate: {length}")
    print(f"unique θ values = {unique_theta_exact}")
    print(f"unique φ values = {unique_phi_exact}")
    print(f"unique ψ values = {unique_psi_exact}")
    print('............................................\n')
    
    #########################################################
    # Apply rotations (one per triplet) and write XYZ
    #########################################################
    # ope file:
    with open(xyz_output_filename, 'w') as file:
        irot = 0
        for theta, phi, psi in triplet_list:   # the rotation is made from triplet
            # make rotation:
            xrot, yrot, zrot = rota(xsv, ysv, zsv, theta, phi, psi, mode)
            
            # write molecule
            file.write(f"{nu}\n")
            file.write(f"{header}\n")
            for iu in range(nu):
                file.write(f"{res[iu]:<8s}"
                           f"{xrot[iu]:10.4f}"
                           f"{yrot[iu]:10.4f}"
                           f"{zrot[iu]:10.4f}\n")
            irot += 1
    
    print('n rotations generated =', irot)
    print(f"XYZ file saved as: {xyz_output_filename}")
    print("angles saved in: data.dat\n")
    print('............................................\n')
    
    # visualize angles:
    filename = f"angles_three_nrot-{irot}_{mode}.png"
    plot_three_angles(theta_list, phi_list, psi_list, filename)
    
    #########################################################
    # Save angles used for plotting on the sphere
    #########################################################
    outfile = "data.dat"
    with open(outfile, "w") as fa:
        for theta, phi, psi in triplet_list:  # the rotation is made from triplet
            fa.write(f"{theta:.10f}  {phi:.10f}  {psi:.10f}\n")
    print("Triplet angles saved as data.dat")
    



########################################################
########################################################

# ===== Runner =====
args = get_cli_args()
input_filename = args.input_file

nu, header, res, x, y, z = read_xyz_file(input_filename)

if args.mode == "grid":
    ntheta = args.ntheta
    nphi   = args.phi
    npsi   = args.psi
    nrot   = ntheta * nphi * npsi
else:
    ntheta = nphi = npsi = None
    nrot   = args.nrot

xyz_output_filename = f"{os.path.splitext(input_filename)[0]}_nrot-{nrot}_{args.mode}.xyz"

print(f"name of the protein: {os.path.splitext(input_filename)[0]}")
adsorbate_rot(
    xyz_output_filename, nu, header, res, x, y, z,
    ntheta, nphi, npsi, args.mode, nrot
)

# Figura de distribución sobre la esfera, usando la convención correcta por modo
angles_to_sphere_points_general(
    input_filename="data.dat",
    mode=args.mode,
    output_filename=f'angles_to_sphere_nrot-{nrot}_{args.mode}.png'
)


########################################################
# un esquema tipico en grilla es:
# ntheta
# nphi = 2*ntheta
# npsi = 2*ntheta
# 6 x 12 x 12 = 864
########################################################


exit()





