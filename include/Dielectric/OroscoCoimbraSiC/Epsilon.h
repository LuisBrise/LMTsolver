// Epsilon_SiC_OC.h -- SiC, Orosco-Coimbra (OC) model. 
//
// omega=0: unlike BB, SiC's OC formula is genuinely analytic here -- no 0/0,
// no pole (confirmed both analytically -- all 11 oscillators have nonzero
// omega0, and mu keeps JiOC's internal "alpha" finite at hbaromega=0 -- and
// numerically in the Python reconstruction: eps_SiC_OC(0) comes out clean,
// ~9.9737, no NaN). So unlike Epsilon_SiC_BB.h's hardcoded 9.9788 and
// Epsilon_Au_OC.h's placeholder 1e10, this file just evaluates the formula
// itself at the near-zero point -- no special constant needed.
//
// VALIDATED against eps_SiC_OC(E_eV) at 10 points, 0.62-30 eV: agreement to
// ~14 significant figures (same eV constant / caveat as the other two files).

std::string material = "OroscoCoimbraSiC";
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

// (kept even though this OC model doesn't use it, for parity with Epsilon.h)
complex<double> JiBB(double hbaromegap, double pesof, double hbarGamma, double hbaromega0, double sigma, double hbaromega){
    return (1i*pesof*pow(hbaromegap,2.0)*sqrt(Pi)/(2*sqrt(2)*hbara(hbarGamma, hbaromega)*sigma))*(Faddeeva::w((hbara(hbarGamma,hbaromega) - hbaromega0)/(sqrt(2.0)*sigma),0.0)+Faddeeva::w((hbara(hbarGamma,hbaromega) + hbaromega0)/(sqrt(2.0)*sigma),0.0));
}

// S(z): analytic continuation used by the OC (Orosco-Coimbra) oscillator. Built
// only from Faddeeva::w -- no other Faddeeva.hh symbol needed for this part.
complex<double> s_of_z(complex<double> z){
    const double eps_reg = 1e-300;
    complex<double> w = Faddeeva::w(z, 0.0);
    double mod_z2_safe = abs(z*z) + eps_reg;
    complex<double> inner = -conj(z) / mod_z2_safe;
    return 1i*Pi*w + exp(-z*z) * (log(z) + log(inner) - 1i*Pi);
}

// OC oscillator. mu is the regularization parameter (1e-5*eV here, matching the
// project's Python convention); "alpha" below is JiOC's own internal variable,
// unrelated to a nanoparticle's polarizability.
complex<double> JiOC(double hbaromegap, double pesof, double hbarGamma, double hbaromega0, double sigma, double hbaromega, double mu){
    double prefactor = pesof * pow(hbaromegap,2.0) / pow(hbaromega0,2.0);
    double Ji0 = -4.0*sqrt(Pi)*Faddeeva::Dawson(-hbaromega0 / (sqrt(2.0)*sigma));
    double sqrt_wG = sqrt(pow(hbaromega,2.0) + pow(hbarGamma,2.0));
    double re_alpha = sqrt(hbaromega/2.0) * sqrt(sqrt_wG + hbaromega);
    double im_alpha = sqrt(hbaromega/2.0) * sqrt(sqrt_wG - hbaromega) + mu;
    complex<double> alpha(re_alpha, im_alpha);
    complex<double> z1 = ( alpha - hbaromega0) / (sqrt(2.0)*sigma);
    complex<double> z2 = (-alpha - hbaromega0) / (sqrt(2.0)*sigma);
    return prefactor * (s_of_z(z1) + s_of_z(z2)) / Ji0;
}

complex<double> eps(double hbaromega)
{
    auto eps_pos = [](double hw)
    {
        const double wp = 10.0 * eV;
        const double mu = 1.0e-5 * eV;

        return 1.0

            // Phonon (shared with the BB file, identical parameters)
            + JiLorentz(
                wp,
                3.184e-4,
                5.893e-4 * eV,
                9.832e-2 * eV,
                hw)

            // OC #1
            + JiOC(wp, 1.081e-1, 2.324e-2 * eV, 3.925    * eV, 3.494    * eV, hw, mu)
            // OC #2
            + JiOC(wp, 3.651e-2, 2.374e-3 * eV, 4.332    * eV, 3.800    * eV, hw, mu)
            // OC #3
            + JiOC(wp, 8.127e-1, 4.164e-4 * eV, 7.128    * eV, 7.740e-1 * eV, hw, mu)
            // OC #4
            + JiOC(wp, 1.462e-1, 5.964e-1 * eV, 6.619    * eV, 3.936e-2 * eV, hw, mu)
            // OC #5
            + JiOC(wp, 3.678e-1, 1.450    * eV, 9.273    * eV, 1.579e-1 * eV, hw, mu)
            // OC #6
            + JiOC(wp, 7.061e-1, 9.884e-1 * eV, 9.217    * eV, 1.466    * eV, hw, mu)
            // OC #7
            + JiOC(wp, 2.089e-1, 1.013e-2 * eV, 1.781e1  * eV, 2.172    * eV, hw, mu)
            // OC #8
            + JiOC(wp, 5.852e-1, 9.285e-3 * eV, 1.292e1  * eV, 2.543    * eV, hw, mu)
            // OC #9
            + JiOC(wp, 5.733e-1, 4.707    * eV, 2.534e1  * eV, 4.734    * eV, hw, mu)
            // OC #10
            + JiOC(wp, 2.615e-1, 2.163    * eV, 1.113e1  * eV, 1.900    * eV, hw, mu)
            // OC #11
            + JiOC(wp, 4.695e-1, 1.183    * eV, 7.326    * eV, 1.831e-1 * eV, hw, mu);
    };

    if (hbaromega > 0.0)
    {
        return eps_pos(hbaromega);
    }
    else if (std::abs(hbaromega) < 1e-14)
    {
        // Genuinely analytic here (see header comment): every sqrt(hbaromega/...)
        // in JiOC multiplies a term that itself -> 0, and no oscillator has
        // omega0=0, so there's no 0/0 anywhere -- safe to evaluate exactly at 0.
        return eps_pos(0.0);
    }
    else
    {
        return conj(eps_pos(-hbaromega));
    }
}

// Omega_SiC_OC.dat -- retuned quadrature partition for SiC, Orosco-Coimbra
// (OC) model. Identical partition to Omega_SiC_BB.dat: BB and OC share the
// exact same phonon term (same Frohlich peak at 0.1155 eV), so the same
// boundaries should work for both.
//
//
// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|
// Omega partition ******************************
const int nw1 = 41;
const int nw2 = 501;
const int nw3 = 281;                   // SiC OC
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
