// Epsilon_SiC_BB.h -- SiC, Brendel-Bormann (BB) model. Same structure/convention
// as Epsilon.h (BrendelBormannAu), just no Drude term (SiC isn't a metal):
// a shared phonon Lorentz term (identical in the OC file too) + 11 BB oscillators,
// Tabla 3.6 of Daniel Espinosa Gonzalez's thesis, transcribed verbatim from his
// Transformada_SiC.cpp (that C++ file disagrees with the thesis TABLE in 3 small
// places.
//
// omega=0: unlike Au, there's no Drude pole, but BB's own hbara() has a genuine
// 0/0 there (confirmed: eps_SiC_BB(0) is NaN in the Python reconstruction too).
// The true static value is finite (~9.979 eV, matching real SiC literature); used
// that (not Au's 1e10, which would be physically wrong here) at the near-zero
// point, same style as Transformada_SiC.cpp's own "return 10.0" (that file used
// a round 10.0; this uses the more precise value, since it's cheap to just say).
//
// VALIDATED against eps_SiC_BB(E_eV) at 10 points, 0.62-30 eV: agreement to ~14
// significant figures (same eV constant / caveat as Epsilon_Au_OC.h).

std::string material = "BrendelBormannSiC";
//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************
#include "Faddeeva.hh"

using namespace std;

//---------- CONSTANTES ----------

//double Pi = boost::math::constants::pi<double>();
double hbar         = 4.135667731/(2*Pi) ;

extern complex<double> Faddeeva::w(complex<double> z, double relerr);

complex<double> hbara(double hbarGamma, double hbaromega){
    return hbaromega*sqrt(1 + sqrt(1.0 + pow(hbarGamma/hbaromega,2.0)))/sqrt(2.0) + 1i*hbaromega*sqrt(-1.0 + sqrt(1.0 + pow(hbarGamma/hbaromega,2.0)))/sqrt(2.0);
}

complex<double> JiLorentz(double hbaromegap, double pesof, double hbarGamma, double hbaromega0, double hbaromega){
    return pow(hbaromegap,2.0)*pesof/(pow(hbaromega0,2.0) - pow(hbaromega,2.0) -1i*hbaromega*hbarGamma);
}

complex<double> JiBB(double hbaromegap, double pesof, double hbarGamma, double hbaromega0, double sigma, double hbaromega){
    return (1i*pesof*pow(hbaromegap,2.0)*sqrt(Pi)/(2*sqrt(2)*hbara(hbarGamma, hbaromega)*sigma))*(Faddeeva::w((hbara(hbarGamma,hbaromega) - hbaromega0)/(sqrt(2.0)*sigma),0.0)+Faddeeva::w((hbara(hbarGamma,hbaromega) + hbaromega0)/(sqrt(2.0)*sigma),0.0));
}

complex<double> eps(double hbaromega)
{
    auto eps_pos = [](double hw)
    {
        const double wp = 10.0 * eV;

        return 1.0

            // Phonon (shared with the OC file, identical parameters)
            + JiLorentz(
                wp,
                3.184e-4,
                5.893e-4 * eV,
                9.832e-2 * eV,
                hw)

            // BB #1
            + JiBB(wp, 3.941e-2, 2.042e-1 * eV, 6.609    * eV, 1.077e-1 * eV, hw)
            // BB #2
            + JiBB(wp, 1.386e-2, 1.426e-4 * eV, 7.965    * eV, 1.885e-1 * eV, hw)
            // BB #3
            + JiBB(wp, 4.353e-2, 4.687    * eV, 3.500    * eV, 2.700e-1 * eV, hw)
            // BB #4
            + JiBB(wp, 1.042e-1, 7.812e-1 * eV, 7.487    * eV, 2.482e-2 * eV, hw)
            // BB #5
            + JiBB(wp, 1.121,    3.305e-1 * eV, 7.101    * eV, 6.155e-1 * eV, hw)
            // BB #6
            + JiBB(wp, 7.599e-2, 9.459e-5 * eV, 1.385e1  * eV, 1.399    * eV, hw)
            // BB #7
            + JiBB(wp, 1.553,    5.689    * eV, 1.032e1  * eV, 1.088    * eV, hw)
            // BB #8
            + JiBB(wp, 8.240e-1, 2.153    * eV, 9.035    * eV, 4.527e-1 * eV, hw)
            // BB #9
            + JiBB(wp, 3.624e-1, 5.376    * eV, 1.691e1  * eV, 1.082    * eV, hw)
            // BB #10
            + JiBB(wp, 3.258e-1, 3.766    * eV, 3.081e1  * eV, 3.973    * eV, hw)
            // BB #11
            + JiBB(wp, 1.556e-2, 3.806e-5 * eV, 9.416    * eV, 2.242e-1 * eV, hw);
    };

    if (hbaromega > 0.0)
    {
        return eps_pos(hbaromega);
    }
    else if (std::abs(hbaromega) < 1e-14)
    {
        // BB's hbara() is 0/0 exactly at hbaromega=0 (genuine, see header comment
        // and sic_causality_pipeline.py); this is the finite physical static
        // limit, not evaluated from the formula.
        return complex<double>(9.9788, 0.0);
    }
    else
    {
        return conj(eps_pos(-hbaromega));
    }
}

// Omega_SiC_BB.dat -- retuned quadrature partition for SiC, Brendel-Bormann
// (BB) model. Same 4-panel Gauss-Kronrod scheme as your WernerAu block (same
// nw/iw formulas, same 5 w's).
//
//
// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|
// Omega partition ******************************
const int nw1 = 41;
const int nw2 = 501;
const int nw3 = 281;                   // SiC BB
const int nw4 = 31;
const int NN = nw1 + nw2 + nw3 + nw4;


double w1 = 0.;
double w2 = 0.003;   // 0.08 eV -- bottom of the phonon window
double w3 = 0.0055;    // 0.14 eV -- top of the phonon window (unchanged from before)
double w4 = 2.0;      // 54.4228 eV
double w5 = 50.;   // wcut (a.u.) -- kept at your own value; still yours to set per b,v.

int iw1 = 2*nw1 + 1;
int iw2 = 2*(nw1 + nw2) + 2;
int iw3 = 2*(nw1 + nw2 +nw3) + 3;
int iw4 = 2*NN + 4;
