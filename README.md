# rotate_adsorbate.py

A Python tool to generate rotated configurations of a molecule or protein using Euler angles (ZYZ convention).  
It supports both grid-based and random sampling and saves the rotated coordinates along with the corresponding angles.  
This tool is useful for molecular simulations, protein orientation studies, and sampling in 3D rotational space.

---

## Requirements

- Python 3.8+
- Dependencies: `numpy`, `scipy`, `matplotlib`

---

## Results

- **Modes**:  
  - `grid` → midpoints (open intervals)  
  - `grid_2` → endpoints included  
  - `random` → Haar-like random sampling  

- **Outputs**:  
  - XYZ trajectory with rotated coordinates  
  - `data.dat` file containing the angles (θ, φ, ψ)  

- **Figures**:  
  - Angle distributions  
  - Sampling preview on the unit sphere  

---

## Uniform sampling of 3D rotations

The rotation of a protein can be represented as the orientation of a unit vector $\hat r$ on the unit sphere ($R = 1$).  
In spherical coordinates:

$$
\hat r = (\sin\theta\cos\phi,\; \sin\theta\sin\phi,\; \cos\theta),
\quad 0 \leq \theta \leq \pi,\; 0 \leq \phi < 2\pi.
$$

The corresponding surface element is:

$$
dS = \sin\theta \, d\theta \, d\phi.
$$

Integrating over the entire sphere:

$$
S = \int_{0}^{\pi} \int_{0}^{2\pi} \sin\theta \, d\theta \, d\phi.
$$

With substitution $u = \cos\theta$:

$$
S = \int_{-1}^{1} \int_{0}^{2\pi} du \, d\phi.
$$

This formulation shows that **uniform sampling in $\cos\theta$ and $\phi$ is required** to generate an unbiased distribution of orientations.  
To account for rotations around $\hat r$, a third angle $\psi$ is introduced with range $0 \leq \psi < 2\pi$.

Thus, a random orientation in 3D is obtained by sampling:

$$
\begin{aligned}
\cos \theta &\in [-1, 1], \\
\phi &\in [0, 2\pi), \\
\psi &\in [0, 2\pi).
\end{aligned}
$$

---

## Quick start

Run the script:

python rotate_adsorbate.py

The script will ask for:
- protein name: e.g. protein.xyz 
- mode: e.g. random or grid
- nrot for random or (ntheta, nphi, npsi) for grid


## Examples results

The script will generate one of the following results in Figure 1.
Panel A shows the distribution of orientations using random sampling with $n_{\text{rot}} = 1000$.
Each dot represents a direction on the unit sphere, defined by the polar angle $\theta$ and the azimuthal angle $\phi$.
It should be noted that $\psi$ angle is distributed similar to $\psi$.
The uniform spread of points across the sphere demonstrates isotropic sampling of orientations.


| Random                                                   | Grid                                                 |
| -------------------------------------------------------- | ---------------------------------------------------- |
| ![Random](figures/angles_to_sphere_nrot-1000_random.png) | ![Grid](figures/angles_to_sphere_nrot-4000_grid.png) |

On the other hand, panel B shows the distribution using grid sampling with $\n_{\theta}$ = 20, $\n_{\phi}$ = 10, and $\n_{\psi}$ = 10.
The figure shows a uniform distribution of the angles on the surface.


On the other hand, the script plot the angles in 2D plot, in order to visualize the extension and distribution of the three applied angles.
The example is shows in grid sampling with with $\n_{\theta}$ = 20, $\n_{\phi}$ = 10, and $\n_{\psi}$ = 10.

![Distribución de ángulos](figures/angles_three_grid.png)

It shoud be noted here, that the a random samplign cover a wide range of angles than grid, however a grid sampling offer a uniform sampling,

## Optional

Additional, from rotate protein file, employing the script 'plot_histo_v6.py', we can calculated  $\n_{\theta}$,  $\n_{\phi}$, and  $\n_{\psi}$ for each generated configuration.
The following figure show the distribution of the angles values with a histogram representation.

![Distribución de ángulos](figures/histos_3x2.png)

This reprentation help to decide the enought amount of notations for the system under investigation.
Ideally, the distribution must be aorund of 0.5 for each configuration.