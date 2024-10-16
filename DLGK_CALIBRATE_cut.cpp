//*****************************************************************************************************************
//
//   DLGK_CALIBRATE.cpp: NOT TESTED
//
//    This program calculates the spectral linear momentum dp/dw transferred to a 
//   NP due a fast swfit electron. It calculate the closed surface integral of
//   the Maxwell stress tensor via Gauss Konrod quadrature, the real part of               
//   the the dpdw void function represents the momentum, and the imaginary part 
//   the estimated error. The integral in frequencies is calculated via Gauss - Kronrod 
//   adaptative algorithm, for low frecuencies, and for high frecuencies an Exp Sinh
//   automatic algorithm is used. Both included in the BOOST libraries.
//   
//
//          By Jesús Castrejon, jcastrejon@ciencias.unam.mx (25/10/2018)
//
//     THIS PROGRAM IS USED FOR CALIBRATING MORE EFICIENT QUADRATURES AND OMEGA INTERVALS!!
//
//*****************************************************************************************************************
//*****************************************************************************************************************

#include <boost/math/quadrature/exp_sinh.hpp>

// Spherical Bessel (Hankel) funcions with *complex argument* (very important) for the scatterred fields. 
// Only found in SLATEC project (FORTRAN 77).

#include <complex_bessel_bits/sph_besselFunctions.h>  

// Source:       https://github.com/valandil/complex_bessel                                                                          

using namespace sp_bessel;
using namespace std;                  

typedef complex<double> dcomplex;          // Defines complex number with double entries.


#include "IN16.h"
#include "EpsilonDrudeAl.h"
#include "AMT_CALIBRATE_cut.h"

//*********************************************************************************************************
//********************************************** MAIN PROGRAM *********************************************
int main(void){

double vv = 0.5;
double a = 5.*nm;                        // impact parameter b, NP ratio, Integration ratio r.
double r = a + 0.05*nm;
double b = a + 0.5*nm;

double lf = 260.;  // 40, v0.3 b1.5;  65, v0.5 b1.5;  100, v0.3 b1.5;  160, v0.9 b1.5; 

auto ff = [a,r,b,vv](double x) { return dl(a,r,b,vv,x) ; };


cout << "Calibrating Al NP" << endl;        // NEEDED TO SET BY HAND !!!!
cout << "v = " << vv << endl;
cout << "b = " << b/nm << endl;
cout << "Lmax = "<< Lmax << endl;
cout << "wc "<< lf << " Har" << endl;
cout << endl;
cout.precision(17);

double l1 = 0.;
double l2 = .2;
double l3 = .4;
double l4 = 2.;
double l5 = std::numeric_limits<double>::infinity();
double termination = sqrt(std::numeric_limits<double>::epsilon());
double error;
double L1;


double Q  = boost::math::quadrature::gauss_kronrod<double, 101>::integrate(ff,  l1, l2, 5, 5.e-14, &error);
double Q1 = boost::math::quadrature::gauss_kronrod<double, 201>::integrate(ff, l2, l3, 5, 5.e-14, &error);
double Q2 = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(ff,  l3, l4, 5, 5.e-14, &error);



boost::math::quadrature::exp_sinh<double> integrator;
double Q3 = integrator.integrate(ff, l4, l5, termination, &error, &L1);


double II = Q + Q1 + Q2 + Q3;


cout << "I = " << II <<endl<< endl;





int i, NN = 24;
double qq[NN], QQ, lcut[NN];

char filename[sizeof "Results/Tail/errorTail_a50nm_b50nm_vv0.99c.dat"];
sprintf(filename, "Results/Tail/errorTail_a%.2gnm_b%.2gnm_vv%.2gc.dat",a/(1.*nm), b/(1.*nm),vv);
FILE *fpp = fopen(filename,"w+");
fprintf(fpp,"Percentual error, a: %.2gnm    v: %.2gc   Lmax: %d  \n", a/(1.*nm), vv, b/(1.*nm), Lmax);
fprintf(fpp,"\n");
fprintf(fpp,"         wcut         error(%)\n");


for ( i = 1; i <= NN ; ++i){

lcut[i] = l4 + (lf-l4)*i/NN;

QQ  = boost::math::quadrature::gauss_kronrod<double, 151>::integrate(ff, l4, lcut[i], 0, 0, &error);

QQ = Q + Q1 + Q2 + QQ;
qq[i] = abs(1. - QQ/II) ;

cout << "{ " << lcut[i] <<","<<  qq[i] << " },"<<endl;
fprintf(fpp,"%.17g %.17g \n",lcut[i], qq[i]);

}


}

//**************************************************** END ****************************************************