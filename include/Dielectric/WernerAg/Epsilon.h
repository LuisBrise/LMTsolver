std::string material = "WernerAg";
//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************
dcomplex LL(double w, double w0, double G, double A) {

  return A/(w0*w0 - w*w - 1i*w*G );
}

dcomplex eps(double w){
return 1. + LL(w,  0.0*eV,     0.065*eV,  98.4*eV*eV)
          + LL(w,  1.0*eV,    93.800*eV, 103.9*eV*eV)
          + LL(w,  4.9*eV,     0.700*eV,  23.8*eV*eV)
          + LL(w, 10.3*eV,    22.000*eV, 318.3*eV*eV)
          + LL(w, 13.2*eV,     6.300*eV,  89.6*eV*eV)
          + LL(w, 21.2*eV,     3.300*eV,  63.5*eV*eV)
          + LL(w, 30.2*eV,     3.700*eV,  43.3*eV*eV)
          + LL(w, 43.0*eV,    16.700*eV, 313.3*eV*eV)
          + LL(w, 65.7*eV,    38.800*eV, 519.9*eV*eV);
}

// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|
// Omega partition ******************************
const int nw1 = 121;   
const int nw2 = 101;   
const int nw3 = 31;                   // Werner Au
const int nw4 = 15;  
const int NN = nw1 + nw2 + nw3 + nw4;


double w1 = 0.;
double w2 = .8;
double w3 = 2.;
double w4 = 5.;
double w5 = 50.;   // This last frequency strongly depends of impact parameter and speed!!
                   // In Au, b=1.5nm a=1nm, for v=0.5 appx 35,  v=0.8 appx  50,v=.99 appx 180..

int iw1 = 2*nw1 + 1;
int iw2 = 2*(nw1 + nw2) + 2;
int iw3 = 2*(nw1 + nw2 +nw3) + 3;
int iw4 = 2*NN + 4;
