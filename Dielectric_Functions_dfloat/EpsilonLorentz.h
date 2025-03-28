#include <string.h>
static char *directoryLv = "LyvLorentz";
static char *directoryLb = "LybLorentz";
static char *directoryPv = "PvLorentz";
static char *directoryPb = "PbLorentz";
//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************
dcomplex LL(double w, double w0, double G, double A) {

  return A/(w0*w0 - w*w - 1i*w*G );
}

dcomplex eps(double w){
return 1. + LL(w, 0.147 , 0.055 , 0.060);
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
