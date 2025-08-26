![cccLogo](https://github.com/user-attachments/assets/96086792-39f4-467c-be08-e3cc3e46b8aa "CCClogo")





# Smoothed Particle Hydrodynamics Simulations – Fluid and Astrophysics🌌

# ![CCCsim logo](https://github.com/user-attachments/assets/769bd30a-241c-4b80-9fa0-56b6213546a2 "CCCsim")
This repository contains 2D Smoothed Particle Hydrodynamics (SPH) simulations written in C++, developed as part of an MSc project at **Panjab University**. The project explores SPH's versatility by applying it to two very different domains: classical fluid dynamics and astrophysical gas interactions.

---

## 🧪 Project Overview

The simulations are divided into two parts:

### 1. **Fluid in a Tumbler (Engineering Simulation)**
A basic test case simulating water-like fluid behavior in a 2D container. This was implemented to validate the SPH method and test pressure dynamics and particle stability under boundary conditions.

### 2. **Cloud-Cloud Collision (Astrophysics Application)**
A more advanced simulation where a stationary V-shaped cloud and a moving rod-shaped cloud collide in 2D. The model tracks how shock compression and cloud interaction lead to the formation of dense clumps—potential precursors to star formation.

---

## ⚙️ Features

- Custom SPH solver written in C++
- 2D particle initialization for both fluid container and cloud geometries
- Pressure, gravity, and artificial viscosity terms included
- Uniform initial density and variable velocity fields
- Parallelized using OpenMP
- Frame-by-frame output for visualization

---

## 🔧 Prerequisites

- **C++11 or newer compiler** (e.g., g++, clang++, or MSVC)
- **OpenMP** support for parallel execution
  - For GCC/MinGW: use `-fopenmp`
  - For MSVC: use `/openmp`
  ```bash
  # To compile 
  g++ -o simulation ccc.cpp -fopenmp -O3 -march=native -ffast-math

- **Make** (optional)
- **Python 3** with `matplotlib` and `numpy` (for visualization)
  ```bash
  pip install matplotlib numpy
---
## 📁 File Structure

```text
├── CCC
    ├── ccc.cpp                    # Main SPH simulation code (cloud-cloud collision)
    ├── ccc_plotter.py             # Python code for visualisation of ccc (to be kept in same directory as ccc.cpp)

├── FLUID_TUMBLER
    ├── fluid_tumbler.cpp          # Optional fluid containment simulation 
    ├── fluid_tumbler_plotter.py   # Python code for visualisation of fluid_tumbler (to be kept in same directory as fluid_tumbler.cpp)

├── visual                         # GIFs of density maps at various time steps are provided in "visual" folder
├── resources                      # various academic resources used for building this project are provided in "resources" folder
