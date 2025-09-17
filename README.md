# rotate_adsorbate.py

Genera configuraciones rotadas de una molécula/proteína usando ángulos de Euler (convención **ZYZ**).  
Soporta muestreos por grilla o aleatorios y guarda tanto las coordenadas como los ángulos usados.

---

## Requisitos

- Python 3.8+
- `numpy`, `scipy`, `matplotlib`

Instalación rápida:
```bash
pip install numpy scipy matplotlib


### Results

- **Modes**: `grid` (midpoints), `grid_2` (endpoints), `random` (Haar-like)
- **Outputs**: XYZ trajectory with rotated coordinates + `data.dat` (θ, φ, ψ)
- **Figures**: angle distributions and sampling preview on the unit sphere


## Quick start

python rotate_adsorbate.py
# Name of the file XYZ (ej: protein.xyz): protein.xyz
# Mode of rotation [grid / random]: random
# Number of divisions in θ: 1000


The following figure shows the angle theta and phi that describe a vector on the surface of an sphere.

![Puntos en la esfera](figures/sphere_points_grid.png)

As we can see, distribution is unifor on the sphere.

Additionaly, the follogin figure shows the range and the distribution of the three angles sampled.
![Distribución de ángulos](figures/angles_three_grid.png)

