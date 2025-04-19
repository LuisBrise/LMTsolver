// include/CoreIncludes.h
#ifndef CORE_INCLUDES_H
#define CORE_INCLUDES_H

// Standard C++ libraries
#include <CLI/CLI.hpp>
#include <filesystem>
#include <string>
#include <iomanip>
#include <sstream>
#include <iostream>    // I/O operations
#include <fstream>     // File operations
#include <cmath>       // Math functions
#include <chrono>      // Timing
#include <complex>     // Complex numbers
#include <vector>      // STL containers (if used)
#include <memory>      // Smart pointers (if used)
#include <quadmath.h>  // To include quadruple precission, requires to compile with -lquadmath

// Spherical Bessel (Hankel) funcions with *complex argument* (very important) for the scatterred fields. 
// Only found in SLATEC project (FORTRAN 77).
#include <complex_bessel_bits/sph_besselFunctions.h> 

using namespace sp_bessel;
using namespace std;
namespace fsm = std::filesystem; 

#endif