# Implementation Guide

This guide explains what the current repository implements, how the code is organized, how it differs from the upstream/original repository, and how to interpret the solver outputs and generated plots.

## 1. Repository Overview

The repository has three practical workflows:

1. **Equilibrium phase-pattern solvers** in `PhaseDiagramCalculations`
   - Four standalone C++ executables:
     - `GradientDescent_1l`
     - `GradientDescent_2l`
     - `SimulatedAnnealing_1l`
     - `SimulatedAnnealing_2l`
2. **Equilibrium post-processing and plotting**
   - Python scripts named `surface_pattern_1l.py` and `surface_pattern_2l.py`
3. **Conserved dynamics on a sphere**
   - A FiPy notebook plus a reusable Python rendering helper in `ConservedDynamics`

The `TraitReconstruction` folder is separate from the solver/plotting pipeline. It contains BayesTraits inputs rather than code that is executed automatically by the C++ or Python workflows.

## 2. Current Implementation Architecture

### 2.1 Build System

The repository now builds from the root `CMakeLists.txt`. This file:

- requires C++11
- finds Eigen automatically when possible
- allows a manual Eigen include path through `POLLEN_EIGEN3_INCLUDE_DIR`
- builds all four solver executables from the root
- optionally builds a Wigner backend test

The canonical build entrypoint is in [CMakeLists.txt](./CMakeLists.txt).

### 2.2 Solver Layout

Each solver is still implemented as a large standalone `main.cpp`. The code is not organized as a modern reusable C++ library. Instead:

- each solver folder contains its own Numerical Recipes headers and local helper types
- most of the numerical workflow lives directly inside `main.cpp`
- model parameters are hardcoded in the source, not passed at runtime

That means each executable represents a specific configured experiment, not a general command-line tool.

### 2.3 Shared Wigner 3j Backend

The current repository includes an in-repo Wigner 3j implementation in:

- `PhaseDiagramCalculations/Common/wigner_compat.h`
- `PhaseDiagramCalculations/Common/wigner_compat.cpp`

This compatibility layer provides the same small API subset the solvers already used:

- `wig_table_init(int, int)`
- `wig_temp_init(int)`
- `wig3jj(int, int, int, int, int, int)`
- `wig_temp_free()`
- `wig_table_free()`

Internally, the implementation:

- checks Wigner 3j selection rules early
- caches log-factorials
- evaluates the Racah-style sum in a numerically stable way
- returns `0.0` immediately for invalid quantum-number combinations

This backend is now the shared mathematical bridge that lets the C++ solvers run without `wigxjpf`.

### 2.4 Python Plotting Layer

The equilibrium plotting scripts:

- load the final spherical-harmonic coefficients from `endcms.txt`
- load the first few displayed parameters and final energy summary from `parameters.txt`
- reconstruct the scalar field on the sphere with `scipy.special`
- normalize the resulting field to `[0, 1]`
- render two views of the sphere using Matplotlib

The current scripts also support both SciPy APIs:

- older `scipy.special.sph_harm`
- newer `scipy.special.sph_harm_y`

That compatibility was added because current SciPy releases removed the older import.

### 2.5 Conserved-Dynamics Layer

The conserved-dynamics workflow now uses:

- `FiPyConservedDynamics.ipynb` for the simulation notebook
- `conserved_dynamics_plot.py` for reusable rendering

The rendering helper:

- converts FiPy mesh face data into 3D coordinates
- optionally displaces the radius by the field amplitude
- renders a 3D Matplotlib scatter plot
- can read bundled FiPy dump files directly

This replaces the old `mayavi`-based plotting dependency.

## 3. The Four Solver Families

## 3.1 `GradientDescent_1l`

This executable minimizes the free-energy/Hamiltonian for a **single spherical-harmonic shell** at one fixed angular degree `l = el_not`.

The active parameters are hardcoded near the top of:

- `PhaseDiagramCalculations/GradientDescent_1l/main.cpp`
- `PhaseDiagramCalculations/GradientDescent_1l/nr.h`

In the current checked-in configuration:

- `el_not = 4`
- `tau = -1`
- `lambda3 = 1`
- `lambda4 = 1`

State representation:

- one coefficient line per `m = 0, 1, ..., l`
- each line stores a real part and an imaginary part
- negative-`m` coefficients are not stored explicitly; they are reconstructed by conjugate symmetry

Optimization method:

- conjugate-gradient style minimization using Numerical Recipes routines
- an analytic derivative is implemented for this case

Outputs:

- `startingcms.txt`
- `cms.txt`
- `endcms.txt`
- `parameters.txt`

In practice, `endcms.txt` is the final solution you plot.

## 3.2 `GradientDescent_2l`

This executable minimizes the Hamiltonian for **two adjacent spherical-harmonic shells**, usually `l = el_not` and `l = el_not + 1`.

In the current source:

- `el_not = 4`
- `kappa = 1`
- `tau = -1`
- `lambda3 = 0`
- `lambda4 = 1`
- `l_not = (2 * el_not + 1) / 2`

Interpretation:

- the model allows coupling between two nearby angular modes
- `kappa` penalizes deviations from the preferred mode center
- `l_not` sets the preferred value around which the quadratic penalty is centered

Optimization method:

- gradient-based minimization again
- but this variant uses numerical differentiation because the full analytic derivative is more cumbersome

Outputs:

- `startingcms.txt`
- `endcms.txt`
- `parameters.txt`
- `time.txt`

The code also opens `cms.txt`, but it is not the main product for interpretation.

## 3.3 `SimulatedAnnealing_1l`

This is the one-`l` version of the model, but optimized with a **simulated annealing / annealed simplex** strategy instead of direct gradient descent.

Current parameters:

- `tau = -1`
- `lambda3 = 0`
- `lambda4 = 1`
- `el_not = 4`

Interpretation:

- same state space as the one-`l` gradient solver
- different optimization strategy
- useful when the energy landscape may have many local minima

Outputs:

- `T_H.txt`
- `T_dH.txt`
- `cms.txt`
- `endcms.txt`
- `parameters.txt`
- `time.txt`

The `T_*.txt` files are annealing diagnostics across the cooling schedule.

## 3.4 `SimulatedAnnealing_2l`

This is the two-shell version optimized by simulated annealing.

Current source settings:

- `l0 = 3`
- `tau = -1`
- `lambda3 = 0`
- `lambda4 = 1`
- `l_not = 3.2`

Interpretation:

- unlike `GradientDescent_2l`, this solver is not parameter-matched to the `el_not = 4` two-shell setup
- it is effectively a different configured experiment

Outputs:

- `T_H.txt`
- `cms.txt`
- `endcms.txt`
- `parameters.txt`
- `time.txt`

## 4. How the Solvers Work

At a high level, all four equilibrium solvers follow the same scientific pipeline:

1. Choose a basis of spherical harmonics.
2. Represent the pattern by complex coefficients in that basis.
3. Build interaction tensors from Gaunt coefficients.
4. Use Wigner 3j symbols to evaluate those Gaunt coefficients.
5. Construct a Hamiltonian/free-energy functional with quadratic, cubic, and quartic terms.
6. Minimize that Hamiltonian numerically.
7. Write the resulting coefficients and diagnostics to text files.
8. Reconstruct the pattern visually with Python.

### 4.1 Role of the Gaunt Coefficients

The Gaunt coefficients encode products of spherical harmonics integrated over the sphere. They are what lets the code translate a field theory on a spherical surface into algebra over harmonic coefficients.

That is why the Wigner backend matters so much: the energy terms are built from `wig3jj(...)`-based Gaunt factors.

### 4.2 Hamiltonian Structure

The parameter names reflect the three main contributions:

- **quadratic term**: controlled by `tau`, and in the 2-`l` models also by `kappa` and `l_not`
- **cubic term**: controlled by `lambda3`
- **quartic term**: controlled by `lambda4`

Very roughly:

- `tau` controls whether patterned states are energetically favored over the trivial state
- `lambda3` controls whether symmetry-breaking cubic interactions are allowed and how strongly they contribute
- `lambda4` stabilizes the amplitude growth through quartic self-interaction
- `kappa` and `l_not` tune the cost of moving away from a preferred angular scale in the mixed-`l` models

## 5. What Changed Relative to the Original Repository

This current repository differs from the upstream/original baseline in several practical ways.

### 5.1 Build and Dependency Changes

Original behavior:

- Linux-style per-folder Makefiles were the primary build path
- `wigxjpf` was required for Wigner 3j values

Current behavior:

- a root `CMakeLists.txt` is now the canonical build path
- `wigxjpf` is no longer required
- the repository ships its own Wigner compatibility implementation

Practical effect:

- the code is much easier to build on Windows
- the numerical dependency chain is shorter

### 5.2 Plotting and Notebook Changes

Original behavior:

- the conserved-dynamics notebook used `mayavi`
- plotting assumptions matched older SciPy

Current behavior:

- `mayavi` is no longer required
- conserved-dynamics rendering uses Matplotlib
- equilibrium plotting scripts support both old and new SciPy spherical-harmonic APIs

Practical effect:

- the repository works better in a modern Python environment
- output can be saved directly to `.png`

### 5.3 Test and Validation Additions

Current repository additions include:

- a Wigner compatibility test
- an optional CTest path
- a smoke-test framework for the solver executables

These were not part of the older workflow.

### 5.4 Local Artifacts vs Source Changes

Folders such as:

- `build/`
- `run_gd1l/`
- generated `.png` files

are run/build artifacts, not part of the conceptual scientific implementation. They are useful examples, but they are not the core repository logic.

## 6. How to Run the Current Workflow

### 6.1 Build

Build once from the repo root:

```bash
cmake -S . -B build
cmake --build build --config Release
```

If Eigen is not auto-detected, point CMake at the directory that contains `Eigen/Dense`.

### 6.2 Run a Solver

Run a solver from the directory where you want its output files to land:

```bash
build/Release/gradient_descent_1l.exe
```

or from a separate run directory:

```bash
mkdir run_gd1l
cd run_gd1l
../build/Release/gradient_descent_1l.exe
```

### 6.3 Plot the Result

For a one-`l` result:

```bash
python PhaseDiagramCalculations/GradientDescent_1l/surface_pattern_1l.py endcms.txt parameters.txt --output gd1l.png
```

For a two-`l` result:

```bash
python PhaseDiagramCalculations/GradientDescent_2l/surface_pattern_2l.py endcms.txt parameters.txt --output gd2l.png
```

Use the script that matches the model family:

- one-`l` solvers -> `surface_pattern_1l.py`
- two-`l` solvers -> `surface_pattern_2l.py`

## 7. Understanding the Output Files

## 7.1 `startingcms.txt`

This is the initial coefficient vector before minimization.

Interpretation:

- it records the starting point in coefficient space
- it is useful for debugging or reproducibility
- it is **not** the final scientific answer

## 7.2 `cms.txt`

This is a trajectory-like output, but the exact meaning depends on the solver.

Interpretation:

- in simulated annealing runs, it often tracks the current best coefficients across the annealing schedule
- in gradient-descent runs, it is less central
- it is usually a diagnostic file rather than the main plotting input

## 7.3 `endcms.txt`

This is the main final coefficient file.

Interpretation:

- it stores the converged spherical-harmonic amplitudes
- it is the primary input for the plotting scripts
- scientifically, this file is the most direct representation of the final equilibrium state

For the one-`l` case, the file stores one line for each `m = 0` through `m = l`, with:

- column 1: real part
- column 2: imaginary part

Example from the current run:

- [run_gd1l/endcms.txt](./run_gd1l/endcms.txt)

## 7.4 `parameters.txt`

This file is both:

- a human-readable run summary
- an input file that the plotting scripts parse for labels

It usually contains:

- the key model parameters
- the starting coefficients
- derivative / convergence information
- Hessian eigenvalues
- the final Hamiltonian value

Example:

- [run_gd1l/parameters.txt](./run_gd1l/parameters.txt)

## 7.5 `T_H.txt` and `T_dH.txt`

These appear in annealing workflows.

Interpretation:

- `T_H.txt`: energy/Hamiltonian across temperature steps
- `T_dH.txt`: derivative-related convergence diagnostic across the cooling schedule

They are useful for optimization diagnostics rather than the final pattern itself.

## 7.6 `time.txt`

This records runtime information when the solver writes it.

Interpretation:

- useful for performance tracking
- not scientifically central to the pattern interpretation

## 8. What the PNG Represents

The equilibrium `.png` file is **not** a raw simulation snapshot. It is a derived visualization built from `endcms.txt` and `parameters.txt`.

### 8.1 What the Plotting Script Does

For the one-`l` case:

1. Read the harmonic coefficients from `endcms.txt`
2. Reconstruct the scalar field on a spherical grid using spherical harmonics
3. Normalize the field to `[0, 1]`
4. Color a sphere with that field
5. Show two viewpoints of the same pattern
6. Add the displayed parameter values and final Hamiltonian value as a title

So the PNG is a visual summary of the final harmonic field on the sphere.

### 8.2 What the Colors Mean

The colors are the normalized value of the reconstructed field.

Interpretation:

- red and blue mark different high/low regions of the reconstructed scalar field
- white or pale transition zones show where the field crosses between those regimes
- the exact sign-to-color mapping is determined by Matplotlib’s `seismic` colormap after normalization

This means the colors are best interpreted **relationally**, not as absolute physical units.

More explicitly, the plotting workflow is:

1. reconstruct \( \phi(\theta,\varphi) \) from the harmonic coefficients
2. find the minimum and maximum of that reconstructed field over the sampled sphere
3. rescale the field values to the interval `[0, 1]`
4. map those normalized values onto the `seismic` colormap

So, within a single PNG:

- blue means the field is near the lower end of that run's sampled values
- red means the field is near the upper end of that run's sampled values
- pale or white regions lie near the midpoint of the normalized range

Two cautions are important:

1. The colors are normalized separately for each run.
   - A red region in one PNG is not guaranteed to represent the same absolute field value as a red region in another PNG.
   - Across runs, the most meaningful comparisons are usually symmetry, patch arrangement, and relative domain structure rather than raw color intensity.

2. White does not necessarily mean \( \phi = 0 \).
   - White is the midpoint between the sampled minimum and maximum after normalization.
   - If the field is not symmetric about zero, that midpoint can differ from the true physical zero of the field.

The safest physical interpretation is therefore:

- red and blue domains are high and low regions of the reconstructed scalar field
- the interfaces between them mark intermediate-value transition zones
- the geometry, connectivity, number, and symmetry of the domains are more robust observables than the exact shade intensity

There is also a useful model-symmetry point:

- when `lambda3 = 0`, the model is often close to a \( \phi \to -\phi \) symmetry
- in that case, swapping red and blue can correspond to the same physical pattern up to sign
- when `lambda3 \neq 0`, that sign symmetry is generally broken more strongly, so polarity can matter more

So the most reliable scientific reading is:

- use the PNG to identify extrema, symmetry class, and domain arrangement
- use `endcms.txt` for the quantitative coefficient state
- use `parameters.txt` and the Hessian/derivative diagnostics to judge convergence and local stability
- avoid assuming that red always means one literal material and blue always means another unless the surrounding scientific interpretation explicitly defines \( \phi \) that way

### 8.3 Why There Are Two Spheres

The one-`l` plotting scripts intentionally render two opposite views of the same sphere.

Interpretation:

- left panel: one viewpoint
- right panel: rotated viewpoint

This helps reveal symmetry and patch arrangement that would be hidden in a single view.

### 8.4 Example: `run_gd1l/gd1l.png`

The current example plot is:

- [run_gd1l/gd1l.png](./run_gd1l/gd1l.png)

It shows:

- a one-`l` (`l = 4`) equilibrium pattern
- strong red/blue contrast over the sphere
- a central blue patch in one view
- paired lateral blue regions in the opposite view

The important thing to read from this image is the **symmetry and arrangement of domains**, not the absolute RGB values.

## 9. Parameter Glossary

## 9.1 `el_not`

Meaning:

- the primary spherical-harmonic angular degree used in the one-`l` model
- the lower of the two degrees in the two-`l` gradient model

Interpretation:

- it sets the angular complexity scale of the basis
- larger values allow finer angular structure on the sphere

## 9.2 Second `l` in the two-`l` model

Meaning:

- the adjacent harmonic shell, typically `el_not + 1`

Interpretation:

- allows the final pattern to mix neighboring angular scales

## 9.3 `tau`

Meaning:

- coefficient of the quadratic term in the Hamiltonian

Interpretation:

- controls whether the zero-amplitude state is favored or whether patterned amplitudes are energetically driven to emerge
- changing `tau` can move the model across a transition between unpatterned and patterned regimes

## 9.4 `lambda3`

Meaning:

- coefficient of the cubic interaction term

Interpretation:

- controls the strength of symmetry-breaking three-mode coupling
- when it is zero, the model has no cubic contribution

## 9.5 `lambda4`

Meaning:

- coefficient of the quartic interaction term

Interpretation:

- stabilizes amplitude growth
- without a stabilizing quartic term, the energy would generally not be bounded in a physically useful way

## 9.6 `kappa`

Meaning:

- quadratic penalty strength in the mixed-`l` model

Interpretation:

- penalizes deviations from the preferred angular mode center
- larger `kappa` forces the solution to stay closer to the preferred scale

## 9.7 `l_not`

Meaning:

- preferred angular mode center in the two-`l` model

Interpretation:

- determines where the quadratic penalty is minimized
- when it lies between two adjacent `l` values, the model can favor mixed-shell behavior

## 9.8 Lowest / Minimum Hamiltonian Value

Meaning:

- the final energy of the converged state

Interpretation:

- lower values are more favorable **within a fixed model and parameter setup**
- you should only compare energies directly when the parameterization and basis are comparable

## 9.9 Hessian Eigenvalues

Meaning:

- local-curvature information of the energy landscape near the converged point

Interpretation:

- positive eigenvalues suggest local stability directions
- near-zero or negative values can indicate flat directions, marginality, or instability

## 10. How to Interpret a Result Scientifically

When you finish a run, the most useful interpretation workflow is:

1. Read the parameter values in `parameters.txt`.
2. Check the final Hamiltonian value.
3. Check whether the reported derivative is small enough to indicate convergence.
4. Look at the Hessian eigenvalues for local stability.
5. Plot `endcms.txt`.
6. Interpret the symmetry and domain structure in the PNG.

Important caution:

- the PNG is a visualization
- the coefficients in `endcms.txt` are the actual quantitative state
- the Hamiltonian and Hessian diagnostics tell you whether that visual state is likely to be a meaningful local minimum

## 11. Example Walkthrough: `run_gd1l`

The current repository includes a concrete example run in:

- [run_gd1l/endcms.txt](./run_gd1l/endcms.txt)
- [run_gd1l/parameters.txt](./run_gd1l/parameters.txt)
- [run_gd1l/gd1l.png](./run_gd1l/gd1l.png)

From `parameters.txt`, the run used:

- `el_not = 4`
- `tau = -1`
- `lambda3 = 1`
- `lambda4 = 1`

The file also reports:

- a starting energy of about `-0.900222`
- `10` iterations
- a very small final derivative
- a final Hamiltonian of about `-15.4905`

Interpretation:

- the optimizer found a state much lower in energy than the starting point
- the small final derivative indicates successful convergence
- the pattern in `gd1l.png` is therefore a plausible converged local minimum for that parameter setting

The `endcms.txt` file contains five stored coefficients because this is the one-`l`, `l = 4` case:

- `m = 0`
- `m = 1`
- `m = 2`
- `m = 3`
- `m = 4`

The plotting script reconstructs the negative-`m` contributions automatically using symmetry, then renders the full spherical field.

## 12. Best Mental Model for the Repository

The simplest accurate way to think about the current repository is:

- the C++ code solves for equilibrium spherical-harmonic coefficients
- the text outputs record those coefficients plus diagnostics
- the Python scripts turn the final coefficient set into an interpretable image
- the conserved-dynamics notebook is a separate, time-dependent workflow
- the current repository modernizes build and plotting compatibility, but the underlying scientific model structure remains the same
