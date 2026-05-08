# CIPHER

**CIPHER (Cauchy-based Inverse Phase Hybrid Error Reduction)** is a fully computational framework to decouple **refractive index (RI)** and **thickness** from a single reconstructed quantitative phase image (QPI), including digital holographic microscopy (DHM) phase maps.

In single-wavelength QPI, the measured phase at each pixel is proportional to the optical path length (OPL) through the specimen:

```math
\phi(x,y)=\frac{2\pi}{\lambda}\,(n(x,y,\lambda)-n_m)\,t(x,y)
```

where $\lambda$ is the illumination wavelength, $n_m$ is the surrounding-medium RI, $n$ is the sample RI, and $t$ is the physical thickness. Because $\phi$ depends on the **product** $(n-n_m)t$, RI and thickness are intrinsically coupled: many $(n,t)$ pairs can produce the same phase value, making direct inversion ill-posed without additional constraints.

CIPHER reduces this ambiguity by incorporating a physically constrained **dispersion model** for the sample RI. In its current implementation, the RI is modeled using a two-parameter Cauchy form:

```math
n(\lambda)=n_1+\frac{n_2}{\lambda^2}
```

which yields the forward phase model:

```math
\phi_{\text{model}}(x,y;n_1,n_2,t)=\frac{2\pi}{\lambda}\left(n_1+\frac{n_2}{\lambda^2}-n_m\right)t
```

Given a reconstructed phase value $\phi_{\text{rec}}$ at a pixel, CIPHER estimates $(n_1,n_2,t)$ by minimizing a phase-mismatch cost (e.g., least squares):

```math
J=(\phi_{\text{rec}}-\phi_{\text{model}})^2
```

subject to physically meaningful bounds on $n_1$, $n_2$, and $t$ and user-defined discretization steps.

## Deterministic vs Optimized CIPHER

### 1) Deterministic CIPHER (exhaustive full-grid)
The deterministic implementation evaluates $J$ on a discrete 3-D lattice spanning $(n_1,n_2,t)$ and selects the minimizer. This guarantees the global minimizer **on the discretized domain**, but the computational cost grows quickly with grid density.

### 2) Optimized CIPHER (FAST full-grid; optional coarse-to-fine refinement)
The optimized implementation accelerates the search by avoiding an explicit loop over thickness. For each $(n_1,n_2)$ grid point, it computes:

```math
\text{denom}=\frac{2\pi}{\lambda}\left(n_1+\frac{n_2}{\lambda^2}-n_m\right),\qquad
t^\*=\frac{\phi_{\text{rec}}}{\text{denom}}
```

then **clamps** $t^\*$ to $[t_{\min},t_{\max}]$ and **quantizes** it to the nearest allowed thickness grid value. This yields the same choice as a brute-force minimization over a *uniform* $t$-grid while reducing runtime dramatically. The method can be embedded in a **coarse-to-fine** strategy: start with coarse steps to locate the solution region, then progressively shrink step sizes and refine within a local window in parameter space.

## Reproducibility and deterministic tie-breaking
In discrete searches, multiple parameter triplets may yield nearly identical minima (within floating-point tolerance), especially in synthetic validations and across CPU/GPU implementations. CIPHER therefore includes an explicit **deterministic tie-breaking rule** (e.g., lexicographic selection of grid indices within a tolerance band) to ensure that the same minimizer is selected reproducibly across runs and platforms.

## Practical notes
- The method is applied **per pixel** (no spatial regularization is required), allowing sharp spatial variations of $n$ and $t$.
- GPU acceleration is natural because cost evaluation over parameter grids is highly parallel.
- Grid step sizes and bounds control the tradeoff between runtime and parameter resolution.

---

This repository contains a minimal, reproducible MATLAB package to replicate **optimized CIPHER** decoupling results (refractive index and thickness) from a **single reconstructed phase map** of a calibrated star phase target.

It includes:
- A MATLAB script that runs the **Optimized CIPHER** implementation on the star target phase map.
- The corresponding `.mat` file containing the star target phase input data.

The goal is to enable readers to reproduce the reported thickness and refractive-index maps (and associated figures/metrics) without requiring access to raw holograms or lab hardware.

---

## Contents

- `run_CIPHER_optimized_star.m`  
  Main script to reproduce results for the star phase target using the Optimized CIPHER implementation.

- `star_input_data.mat`  
  Input dataset used by `run_CIPHER_optimized_star.m`. This file contains the reconstructed phase map (in radians).

---

## Requirements

**MATLAB** (R2019b or newer recommended)

Toolboxes (recommended):
- **Parallel Computing Toolbox** (for GPU acceleration via `gpuArray`, if enabled)
- **Image Processing Toolbox** 

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
2. Defines physical parameters (e.g., wavelength $\lambda$ and surrounding refractive index $n_m$).
3. Sets search ranges and grid steps for the CIPHER parameters $(n_1, n_2, t)$.
4. Runs the **Optimized CIPHER** inversion:
   - Uses a discrete search over $(n_1, n_2)$
   - Computes the best thickness $t$ by fast selection on the discretized grid (FAST full-grid logic)
   - Uses deterministic tie-breaking to ensure reproducibility when multiple minimizers exist.
5. Displays and/or saves:
   - Thickness map $t(x,y)$  
   - Refractive index map $n(x,y)$ at the illumination wavelength  
   - Optional intermediate maps such as $n_1(x,y)$, $n_2(x,y)$, residual/cost maps, and filtered versions.

---

## Inputs / Outputs

### Inputs
- `star_input_data.mat` (loaded by the script)

Typical variables inside the script include:
- `phi_map` (reconstructed phase map in radians loaded from .mat file)
- `lambda_um` or `lambda` (wavelength)
- `nm` (surrounding medium refractive index)
- Optional masks or reference maps (depending on your implementation)

### Outputs
By default, the script places results in the MATLAB workspace and produces figures.

If your script saves outputs (recommended), we suggest saving:
- `t_map`, `n_map`, `n1_map`, `n2_map`, `res2_map`
- Any filtered variants: `t_map_fil`, `n_map_fil`
- A results `.mat` file (e.g., `outputs/star_results.mat`)

---

## Reproducibility Notes

- **GPU vs CPU:** Results should match numerically (within floating tolerance), but exact tie cases can produce different selected indices if deterministic tie-breaking is not enforced. This script includes deterministic tie-breaking so that the output matches a brute-force “first-hit” strategy consistently.
- **Grid step sizes matter:** $\Delta n_1$, $\Delta n_2$, and $\Delta t$ directly affect both accuracy and runtime. The script uses the same settings reported in the manuscript to reproduce the figures.

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

If you use this code or data, please cite our manuscript (under revision in Elsevier's  Measurement Journal):

C. Joseph, C. Trujillo, and A. Doblas, “Metrological decoupling of refractive index and thickness in quantitative phase imaging using dispersion-based metrological modeling,” submitted for publication in Measurement (Elsevier), 2026.

**Title of the Manuscript:**  
*Metrological Decoupling of Refractive Index and Thickness in Quantitative Phase Imaging using Dispersion-Based Metrological Modeling*

**Authors:**  
Clivens Joseph¹, Carlos Trujillo², Ana Doblas¹

**Affiliations:**  
¹ Department of Electrical and Computer Engineering, University of Massachusetts Dartmouth, Dartmouth, MA, United States  
² Optics and Photonics Laboratory, School of Applied Science and Engineering, Universidad EAFIT, Medellín, Colombia

**Author e-mails:**  
Ana Doblas — adoblas@umassd.edu
Carlos Trujillo — catrujilla@eafit.edu.co
