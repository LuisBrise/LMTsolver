#include <string.h>
static char *directoryLv = "LyvAu";
static char *directoryLb = "LybAu";
static char *directoryPv = "PvAu";
static char *directoryPb = "PbAu";
//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************
dcomplex LL(double w, double w0, double G, double A) {

  return A/(w0*w0 - w*w - 1i*w*G );
}

dcomplex eps(double w){
return 1. + LL(w, 0.                , 0.007349864496171 , 0.152742986686889)
          + LL(w, 0.146997289923419 , 0.055123983721282 , 0.060232866544963)
          + LL(w, 0.268270054110239 , 0.12127276418682  , 0.074008096113541)
          + LL(w, 0.47039132775494  , 0.433642005274085 , 0.249709798748062)
          + LL(w, 0.694562194888153 , 2.60920189614068  , 0.983308298910026)
          + LL(w, 0.731311517369008 , 0.106573035194478 , 0.088728684574082)
          + LL(w, 1.0620554196967   , 0.143322357675333 , 0.067525635140093)
          + LL(w, 1.42219878000908  , 0.477741192251111 , 0.100883298899298)
          + LL(w, 2.36298143551895  , 1.90728983675636  , 0.734678910324206);
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
