// Epsilon_Si_Lorentz.h -- Si, Drude-Lorentz model. 

std::string material = "LorentzSi";
//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************

using namespace std;

inline dcomplex LorentzPoleSi(double w,
                              double wp, double wpp,
                              double sp, double spp)
{
    const dcomplex Omega(wp, wpp);
    const dcomplex sigma(sp, spp);

    // Explicit conjugates, so this does not depend on std::conj for dcomplex.
    const dcomplex sigma_conj(sp, -spp);
    const dcomplex Omega_conj(wp, -wpp);

    const dcomplex I(0.0, 1.0);

    // Paper Eq. (2):
    // i*sigma/(w - Omega) + i*sigma*/(w + Omega*)
    return I*sigma/(w - Omega) + I*sigma_conj/(w + Omega_conj);
}

inline dcomplex eps(double w)
{
    // Si parameters from Table 1 of the paper:
    // epsilon_inf = 0.81568
    // Omega_k = Omega'_k + i Omega''_k
    // sigma_k = sigma'_k + i sigma''_k
    return 0.81568
        + LorentzPoleSi(w, 3.37361*eV, -0.11402*eV,  1.6934*eV,  2.0840*eV)
        + LorentzPoleSi(w, 3.65190*eV, -0.52378*eV,  5.2573*eV,  8.0106*eV)
        + LorentzPoleSi(w, 4.28770*eV, -0.21116*eV, -1.7164*eV,  5.9939*eV)
        + LorentzPoleSi(w, 5.31880*eV, -0.18434*eV, -0.00528*eV, 0.32911*eV)
        + LorentzPoleSi(w, 5.50640*eV, -1.78920*eV, -3.8438*eV,  6.9298*eV);
}

// Omega_SiC_BB.dat -- retuned quadrature partition for SiC, Brendel-Bormann
// (BB) model. Same 4-panel Gauss-Kronrod scheme as your WernerAu block (same
// nw/iw formulas, same 5 w's).
//
// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|
// Omega partition ******************************
const int nw1 = 41;
const int nw2 = 501;
const int nw3 = 281;                   // Si Lorentz
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
