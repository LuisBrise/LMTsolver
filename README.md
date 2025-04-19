# AMTsolver - Angular Momentum Transfer Solver

## Overview

AMTsolver is a C++ program that calculates the spectral angular momentum (dL/dω) transferred to a nanoparticle (NP) due to interaction with a fast-moving electron. The calculation is performed by analytically solving the closed surface integral, resulting in a double multipolar sum expression. The program focuses exclusively on the external field contribution and calculates the total angular momentum transfer (ΔL) for different values of the impact parameter b.

## Features

- Computes frequency integrals using Gauss-Kronrod quadrature with adaptive partition
- Saves spectral results to `dldw_*.dat` files
- Saves integrated angular momentum results to `DL_*.dat` files
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
git clone https://github.com/yourusername/AMTsolver.git
cd AMTsolver
mkdir build

## Basic compilation
g++ -O2 -std=c++17 -Iinclude -c src/main.cpp -o build/main.out -lcomplex_bessel
g++ build/IN51.out build/main.out -o AMTsolverL51 -lcomplex_bessel

## With OpenMP support
g++ -fopenmp -O2 -std=c++17 -Iinclude -c src/main.cpp -o build/main.out -lcomplex_bessel
g++ build/IN51.out build/main.out -o AMTsolverL51 -lcomplex_bessel -fopenmp

## Basic Command (example)
./AMTsolverL51 -a <radius_nm> -b <impact_param_nm> -v <velocity>

## Single Calculation
./AMTsolverL51 -a 50.0 -b 55.0 -v 0.7

## Velocity Scan (20 points)
./AMTsolverL51 -a 50.0 -b 55.0 --vscan -n 20
## Impact parameter scan
./AMTsolverL51 -a 50.0 -v 0.7 --bscan -n 30 -e 1e-5

## Command Line Options

| Option       | Description                                      | Default Value |
|--------------|--------------------------------------------------|---------------|
| `-a`         | Nanoparticle radius (nm)                         | 1.0           |
| `-b`         | Impact parameter (nm)                            | 1.5           |
| `-v`         | Electron velocity (fraction of speed of light)   | 0.7           |
| `--vscan`    | Enable velocity scan mode                        | false         |
| `--bscan`    | Enable impact parameter scan mode                | false         |
| `--contour`  | Enable velocity vs. impact parameter contour     | false         |
| `-e`         | Convergence threshold for multipolar expansion   | 1e-4          |
| `-n`         | Number of points in scans                        | 9             |
| `-m`         | Minimum multipole order (ℓ) for calculations     | 1             |

**Notes:**
- The impact parameter (`-b`) must be greater than the NP radius (`-a`)
- Velocity (`-v`) must be between 0.0 and 0.999 (in units of c)
- When using scan modes (`--vscan` or `--bscan`), the `-n` option controls the number of points

## Output structure
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
            
## For citations
@software{AMTsolver,
  author = {Jorge Luis Briseño-Gómez},
  title = {AMTsolver: Angular Momentum Transfer Calculator},
  year = {2023},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/yourusername/AMTsolver}}
