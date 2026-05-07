# Optimized CIPHER (MATLAB)

This repository contains a minimal, reproducible MATLAB package to replicate optmized **CIPHER** decoupling results (refractive index and thickness) from a **single reconstructed phase map** of a calibrated star phase target.

It includes:
- A MATLAB script that runs the **Optimized CIPHER** implementation on the star target phase map.
- The corresponding `.mat` file containing the star target phase input data.

The goal is to enable readers to reproduce the reported thickness and refractive-index maps (and associated figures/metrics) without requiring access to raw holograms or lab hardware.

---

## Contents

- `run_CIPHER_optimized_star.m`  
  Main script to reproduce results for the star phase target using the Optimized CIPHER implementation.

- `star_input_data.mat`  
  Input dataset used by `run_CIPHER_optimized_star.m`. This file contains the reconstructed phase map (in radians) and any auxiliary parameters required by the script (e.g., wavelength, medium refractive index, etc., depending on how the script is set up).

---

## Requirements

**MATLAB** (R2019b or newer recommended)

Toolboxes (recommended):
- **Parallel Computing Toolbox** (for GPU acceleration via `gpuArray`, if enabled)
- **Image Processing Toolbox** (if the script uses filtering, thresholding, or segmentation)

Hardware:
- GPU acceleration is optional. The code runs on CPU as well (slower).

---

## Quick Start

1) Clone or download this repository:
```bash
git clone https://github.com/catrujilla/CIPHER.git
cd CIPHER
```

2) Open MATLAB and set the repo folder as your working directory:
```matlab
cd('<path_to_repo>')
```

3) Run the script:
```matlab
run('run_CIPHER_optimized_star.m')
```

---

## What the script does

At a high level, `run_CIPHER_optimized_star.m`:

1. Loads the star phase input data from `star_input_data.mat`.
2. Defines physical parameters (e.g., wavelength λ and surrounding refractive index \(n_m\)).
3. Sets search ranges and grid steps for the CIPHER parameters \((n_1, n_2, t)\).
4. Runs the **Optimized CIPHER** inversion:
   - Uses a discrete search over \((n_1, n_2)\)
   - Computes the best thickness \(t\) by fast selection on the discretized grid (FAST full-grid logic)
   - Uses deterministic tie-breaking to ensure reproducibility when multiple minimizers exist.
5. Displays and/or saves:
   - Thickness map \(t(x,y)\)  
   - Refractive index map \(n(x,y)\) at the illumination wavelength  
   - Optional intermediate maps such as \(n_1(x,y)\), \(n_2(x,y)\), residual/cost maps, and filtered versions.

---

## Inputs / Outputs

### Inputs
- `star_input_data.mat` (loaded by the script)

Typical variables inside may include:
- `phi_map` (reconstructed phase map in radians)
- `lambda_um` or `lambda` (wavelength)
- `nm` (surrounding medium refractive index)
- Optional masks or reference maps (depending on your implementation)

### Outputs
By default, the script places results in the MATLAB workspace and produces figures.

If your script saves outputs (recommended), we suggest saving:
- `t_map`, `n_map`, `n1_map`, `n2_map`, `res2_map`
- Any filtered variants: `t_map_fil`, `n_map_fil`
- A results `.mat` file (e.g., `outputs/star_results.mat`)
- Exported figures (e.g., PNG/SVG/PDF)

---

## Reproducibility Notes

- **GPU vs CPU:** Results should match numerically (within floating tolerance), but exact tie cases can produce different selected indices if deterministic tie-breaking is not enforced. This script includes deterministic tie-breaking so that the output matches a brute-force “first-hit” strategy consistently.
- **Grid step sizes matter:** \(\Delta n_1\), \(\Delta n_2\), and \(\Delta t\) directly affect both accuracy and runtime. The script uses the same settings reported in the manuscript to reproduce the figures.

---

## Troubleshooting

**1) `Undefined function 'gpuArray'...`**  
You need the **Parallel Computing Toolbox**, or disable GPU in the script:
```matlab
opts.use_gpu = false;
```

**2) Script runs but values appear “stuck” at bounds**  
This may indicate grid ambiguity or low-information regions. Check:
- segmentation/mask (exclude background)
- bounds/step sizes
- deterministic tie-break tolerance

**3) No figures / no outputs saved**  
Check whether the script is configured to save. If not, add a section like:
```matlab
save('star_results.mat','t_map','n_map','n1_map','n2_map','res2_map','opts','seeds');
```

---

## How to Cite

If you use this code or data, please cite our manuscript (under revision):

**Title of the Manuscript:**  
*Metrological Decoupling of Refractive Index and Thickness in Quantitative Phase Imaging using Dispersion-Based Metrological Modeling*

**Authors:**  
Clivens Joseph¹, Carlos Trujillo², Ana Doblas¹

**Affiliations:**  
¹ Department of Electrical and Computer Engineering, University of Massachusetts Dartmouth, Dartmouth, MA, United States  
² Optics and Photonics Laboratory, School of Applied Science and Engineering, Universidad EAFIT, Medellín, Colombia

**Corresponding Author:**  
Ana Doblas — adoblas@umassd.edu

Carlos Trujillo — catrujilla@eafit.edu.co
