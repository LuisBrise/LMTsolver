#include <string.h>
static char *directoryLv = "LyvDrude";
static char *directoryLb = "LybDrude";
static char *directoryPv = "PvDrude";
static char *directoryPb = "PbDrude";
//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************
//*********************** parameters given by Rakic ************************
double wp     = 0.367493;
double Gamma  = 0.00367493;
dcomplex eps(double w){
return 1. - pow(wp,2.)/(w*(w + 1i*Gamma));
} 

// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|
const int nw1 = 51;   
const int nw2 = 301;   
const int nw3 = 51;                   // Drude Al
const int nw4 = 51;  
const int NN = nw1 + nw2 + nw3 + nw4;

double w1 = 0.;
double w2 = .05;
double w3 = .3;
double w4 = 2.;
double w5 = 40.; // Omega cut for freq. integral

int iw1 = 2*nw1 + 1;
int iw2 = 2*(nw1 + nw2) + 2;
int iw3 = 2*(nw1 + nw2 +nw3) + 3;
int iw4 = 2*NN + 4;
