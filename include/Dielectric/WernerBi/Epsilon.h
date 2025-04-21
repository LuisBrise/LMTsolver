std::string material = "WernerBi";
//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************
dcomplex LL(double w, double w0, double G, double A) {

  return A/(w0*w0 - w*w - 1i*w*G );
}

dcomplex eps(double w){
return 1. + LL(w,  0.0*eV,     0.087*eV, 108.4*eV*eV)
          + LL(w,  3.9*eV,     1.800*eV,  23.6*eV*eV)
          + LL(w,  6.6*eV,     4.500*eV,  50.2*eV*eV)
          + LL(w, 11.3*eV,     6.200*eV,  44.8*eV*eV)
          + LL(w, 24.1*eV,     0.600*eV,   6.6*eV*eV)
          + LL(w, 27.5*eV,     1.900*eV,  13.2*eV*eV)
          + LL(w, 40.2*eV,     9.600*eV,  99.3*eV*eV)
          + LL(w, 55.8*eV,    21.700*eV, 326.8*eV*eV)
          + LL(w, 84.2*eV,    48.400*eV, 373.1*eV*eV);
}

// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|
// Omega partition ******************************
const int nw1 = 121;   
const int nw2 = 181;   
const int nw3 = 101;                   // Werner Au
const int nw4 = 31;  
const int NN = nw1 + nw2 + nw3 + nw4;


double w1 = 0.;
double w2 = .22;
double w3 = 1.;
double w4 = 5.;
double w5 = 100.;   // This last frequency strongly depends of impact parameter and speed!!
                   // In Au, b=1.5nm a=1nm, for v=0.5 appx 35,  v=0.8 appx  50,v=.99 appx 180..

int iw1 = 2*nw1 + 1;
int iw2 = 2*(nw1 + nw2) + 2;
int iw3 = 2*(nw1 + nw2 +nw3) + 3;
int iw4 = 2*NN + 4;
