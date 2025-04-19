// include/CoreTypes.h
#ifndef CORE_TYPES_H
#define CORE_TYPES_H

#include <complex>
//#include "mpreal.h"

//mpfr::mpreal::set_default_prec(256);
// Fundamental constants and types
constexpr double Cspeed = 137.035999139;
constexpr double nm = 100./5.2918;
constexpr double eV = 1.0/27.211386245988;
//using dcomplex = std::complex<double>;  // Modern C++ alternative to typedef
typedef std::complex<double> dcomplex;
//typedef mpfr::mpreal qfloat;
//typedef std::complex<qfloat> qcomplex;

struct SimulationConfig {
    std::string material;
    int Lmax;
    double a;        // NP radius [nm]
    double b;        // Impact parameter [nm]
    double velocity; // [c units]
    const bool isVScan;
    const bool isBScan;
    const bool isBvsVContour;
    double bInit; // Initial value of impact parameter in scan
    double bFin;  // Final value of impact parameter in scan 
    double vvInit;  // Initial value of speed in scan
    double vvFin;   // Final value of speed in scan
    const double errorThreshold;
    const int numOfPoints; // Number of points for b
    int minLmax; // Minimum multipolar order
    std::string timestamp;
    double r;  // Always a + 0.05 (enforced at construction)    

    // Constexpr constructor for compile-time checks
    SimulationConfig(std::string mat, int lmax, 
                     double np_radius, double impact_param, double v,
                     bool vscan, bool bscan, bool bvsv,
                     double bi, double bf, double vi, double vf,
                     double err_thresh, int num_points, int minL,
                     std::string tmstmp)
        : material(std::move(mat)),
          Lmax(lmax),
          a(np_radius),
          b(impact_param),
          velocity(v),
          isVScan(vscan),
          isBScan(bscan),
          isBvsVContour(bvsv),
          bInit(bi),
          bFin(bf),
          vvInit(vi),
          vvFin(vf),
          errorThreshold(err_thresh),
          numOfPoints(num_points),
          minLmax(minL),
          timestamp(tmstmp),
          r(np_radius + 0.05)  // Enforced relationship 
          { 
        if (a <= 0) throw std::invalid_argument("NP radius must be positive");
    }
};  

#endif