#include <string.h>
static char *directoryLv = "LyvAl";
static char *directoryLb = "LybAl";
static char *directoryPv = "PvAl";
static char *directoryPb = "PbAl";
//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************
//*********************** parameters given by Rakic ************************
dfloat wp     = 0.4829593006;
dfloat Gamma  = 0.0072396121;
dcomplex eps_dfloat(dfloat w){
return 1. - pow(wp,2.)/(w*(w + imag_unit*Gamma));
} 

// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|
const int nw1 = 15;   

dfloat w1 = 0.;
dfloat w2 = .1;

// Parameters to calculate Gauss Kronrod weights

int iw1 = 2*nw1 + 1;

/*
********* Use the following parameters for a=5nm, 10nm
const int nw1 = 51;   
const int nw2 = 201;   
const int nw3 = 51;                   // Drude Al
const int nw4 = 51;  
const int NN = nw1 + nw2 + nw3 + nw4;

double w1 = 0.;
double w2 = .2;
double w3 = .4;
double w4 = 2.;
double w5 = 40.; // Omega cut for freq. integral
*/

/*
********* Use the following parameters for a=20 nm, 50 nm
const int nw1 = 51;   
const int nw2 = 251;   
const int nw3 = 51;                   // Drude Al
const int nw4 = 51;  
const int NN = nw1 + nw2 + nw3 + nw4;

double w1 = 0.;
double w2 = .15;
double w3 = .4;
double w4 = 2.;
double w5 = 40.; // Omega cut for freq. integral
*/