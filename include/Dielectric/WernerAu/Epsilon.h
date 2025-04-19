std::string material = "WernerAu";
//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************
dcomplex LL(double w, double w0, double G, double A) {

  return A/(w0*w0 - w*w - 1i*w*G );
}

dcomplex eps(double w){
return 1. + LL( 0.0*eV,     0.2*eV,   113.1*eV*eV)
          + LL( 4.0*eV,     1.5*eV,    44.6*eV*eV)
          + LL( 7.3*eV,     3.3*eV,    54.8*eV*eV)
          + LL(12.8*eV,    11.8*eV,   184.9*eV*eV)
          + LL(18.9*eV,    71.0*eV,   728.1*eV*eV)
          + LL(19.9*eV,     2.9*eV,    65.7*eV*eV)
          + LL(28.9*eV,     3.9*eV,    50.0*eV*eV)
          + LL(38.7*eV,    13.0*eV,    74.7*eV*eV)
          + LL(64.3*eV,    51.9*eV,   544.0*eV*eV)
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
