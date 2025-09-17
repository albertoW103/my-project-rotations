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


- **Modes**: `grid` (midpoints), `grid_2` (endpoints), `random` (Haar-like)
- **Outputs**: XYZ trajectory with rotated coordinates + `data.dat` (θ, φ, ψ)
- **Figures**: angle distributions and sampling preview on the unit sphere


## Resultados


### Puntos sobre la esfera
![Puntos en la esfera](figures/sphere_points_grid.png)


### Distribución de ángulos (θ, φ, ψ)
![Distribución de ángulos](figures/angles_three_grid.png)

## Quick start
python rotate_adsorbate.py
# Name of the file XYZ (ej: protein.xyz): protein.xyz
# Mode of rotation [grid / random]: grid
# Number of divisions in θ: 6
# Number of divisions in φ: 12
# Number of divisions in ψ: 12