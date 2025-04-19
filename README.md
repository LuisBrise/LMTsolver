# LMTsolver - Linear Momentum Transfer Solver

## Overview

LMTsolver is a C++ program that calculates the spectral linear momentum (dP/dω) transferred to a nanoparticle (NP) due to interaction with a fast-moving electron. The calculation is performed by analytically solving the closed surface integral, resulting in a double multipolar sum expression. The program focuses exclusively on the external field contribution and calculates the total angular momentum transfer (ΔP) for different values of the impact parameter b.

## Features

- Computes frequency integrals using Gauss-Kronrod quadrature with adaptive partition
- Saves spectral results to `dpdw_*.dat` files
- Saves integrated angular momentum results to `DP_*.dat` files
- Supports multiple parameter sets via CLI
- Automatic convergence checking with configurable threshold
- Parallel computation support (via OpenMP)

## Prerequisites

- C++17 compatible compiler (g++ recommended)
- Boost libraries (version 1.70+)
- complex_bessel library
- CLI11 (for command line parsing)

## Installation

```bash
git clone https://github.com/yourusername/LMTsolver.git
cd LMTsolver
mkdir build
```

# Standard build
```bash
g++ -Iinclude -c src/IN51.cpp -o build/IN51.o
g++ -Iinclude -Iinclude/Dielectric/DrudeAl -c src/main.cpp -o build/main.o -lcomplex_bessel
g++ build/IN51.o build/main.o -o bin/AMTsolver -lcomplex_bessel
```

# Parallel build (OpenMP)
```bash
g++ -fopenmp -O2 -std=c++17 -Iinclude -c src/main.cpp -o build/main.o -lcomplex_bessel
g++ build/IN51.o build/main.o -o bin/AMTsolver -lcomplex_bessel -fopenmp```
```

## Basic Command (example)
```bash
./LMTsolver -a <radius_nm> -b <impact_param_nm> -v <velocity>
```

## Single Calculation
```bash
./LMTsolver -a 50.0 -b 55.0 -v 0.7
```
## Velocity Scan (20 points)
```bash
./LMTsolver -a 50.0 -b 55.0 --vscan -n 20
```
## Impact parameter scan
```bash
./LMTsolver -a 50.0 -v 0.7 --bscan -n 30 -e 1e-5
```
## Command Line Options

| Option       | Description                                      | Default Value | Constraints                      |
|--------------|--------------------------------------------------|---------------|----------------------------------|
| `-a <float>` | Nanoparticle radius (nm)                         | 1.0           | Must be > 0                      |
| `-b <float>` | Impact parameter (nm)                            | 1.5           | Must be > NP radius (`-a`)       |
| `-v <float>` | Electron velocity (fraction of c)                | 0.7           | 0.0 ≤ v ≤ 0.999                  |
| `--vscan`    | Enable velocity scan mode                        | false         | Requires `-n` for point count    |
| `--bscan`    | Enable impact parameter scan mode                | false         | Requires `-n` for point count    |
| `--contour`  | Enable velocity vs. impact parameter contour     | false         | Combines `--vscan` and `--bscan` |
| `-e <float>` | Convergence threshold                            | 1e-4          | Must be > 0                      |
| `-n <int>`   | Number of points in scans                        | 9             | Must be ≥ 2                      |
| `-m <int>`   | Minimum multipole order (ℓ)                      | 1             | Must be ≥ 1                      |

**Notes:**
- All numerical values must be positive
- Scans generate data files in the output directory structure
- Contour mode generates 2D parameter space analysis

## Output structure
```text
results/
└── material=<material>/
    └── a=<radius>nm/
        └── <timestamp>/
            ├── v_scan_for_b=<impact_param>nm/  # Velocity scans
            │   ├── Lmax=<value>/
            │   │   ├── DLy_<params>.dat
            │   │   └── error_DLy_<params>.dat
            │   └── MultipolarConvergence/
            │       └── MultipolarConvergence_<params>.dat
            ├── b_scan_for_v=<velocity>c/       # Impact param scans
            │   └── ... 
            └── simulation_parameters.txt       # Metadata
```
            
## For citations
@software{AMTsolver,
  author = {Jesús Castrejón-Figeroa and Jorge Luis Briseño-Gómez},
  title = {LMTsolver: Linear Momentum Transfer Calculator},
  year = {2023},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/xuxocast/AMTsolver}}
