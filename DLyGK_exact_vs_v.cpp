//*****************************************************************************************************************
//
//   DLGK_exact_vs_b.cpp:
//
//    This program calculates the spectral linear momentum dl/dw transferred to a 
//   NP due a fast swfit electron. It calculate the closed surface integral 
//   analitically, obteining a double multipolar sum espresion. The spectrum is 
//   written to a file named "dldw*.dat". The program ONLY considers the external field.
//   Also, the integral in frequencies is calculated via Gauss - Kronrod using
//   a partition set to fast convergence. The results are written to a file named "DL*.dat". 
//  
//	flags 
//      g++ -o DLGK_exact_vs_v.out DLGK_exact_vs_v.cpp -lcomplex_bessel -w 
//      g++ -o DLyGK_exact_vs_v.out DLyGK_exact_vs_v.cpp -lcomplex_bessel -w && ./DLyGK_exact_vs_v.out 
//
//          By Jesús Castrejon, jcastrejon@ciencias.unam.mx (25/02/2019)
// Modified by Jorge Luis Briseño, jorgeluisbrisenio@ciencias.unam.mx (21/02/2023)
//
//*****************************************************************************************************************
//*****************************************************************************************************************
// Spherical Bessel (Hankel) funcions with *complex argument* (very important) for the scatterred fields. 
// Only found in SLATEC project (FORTRAN 77).

#include <complex_bessel_bits/sph_besselFunctions.h>  

// Source:       https://github.com/valandil/complex_bessel                                                                          

using namespace sp_bessel;
using namespace std;                  

typedef complex<double> dcomplex;          // Defines complex number with double entries.

char directory[10] = "LyvAl";
#include "EpsilonDrudeAl.h"
#include "AngularMomentumTransfer_l31.h"
//*********************************************************************************************************
//********************************************** MAIN PROGRAM *********************************************
int main(void){

cout.precision(17);
cout << endl;
cout << "Momentum Gauss - Kronrod :" << endl;
cout << endl;
//cout << "Lmax = " << Lmax << endl;
//cout << endl;

/*
******** The selection of parameters for a sweep ********
******** in the speed of the electron. Non used *********
******** parameters are just commented. *****************
*/

/*double a = 5.0*nm; 
double b = a + 0.5*nm;                                               
double r = a + 0.05*nm;
double vv = 0.5;*/

// Radius of NP 10, 20, 50
double a = 50.0*nm;
double b = a + 1.0*nm;                                                
double r = a + 0.05*nm;
double vv = 0.5;

dcomplex DLy[6];
double Ly;

// Start the timer
auto start = std::chrono::high_resolution_clock::now();
for (Lmax = 10; Lmax < 31; Lmax++){
cout << "Lmax = " << Lmax << endl;
cout << endl;
char filename[sizeof "Results/DL_a1nm_v0.99c_b1.5nm_extF2.dat"];
sprintf(filename, "Results/%s/%d/DL_a%.2gnm_b%.2gnm.dat", directory,Lmax,a/(1.*nm), b/(1.*nm));
FILE *fpp = fopen(filename,"w+");
fprintf(fpp,"Total Angular momentum transfered, a: %.2gnm    v: %.2gc   Lmax: %d  \n", a/(1.*nm), vv, b/(1.*nm), Lmax);
fprintf(fpp,"\n");
fprintf(fpp,"         betav         DLEy                  DLHy                  DLEsy                  DLHsy                  DLEsyExt                  DLHsyExt                  DLy\n");

char filenamer[sizeof "Results/DL_a1nm_v0.99c_b1.5nm_error_extF2.dat"];
sprintf(filenamer, "Results/%s/%d/DL_a%.2gnm_b%.2gnm_error.dat", directory,Lmax, a/(1.*nm), b/(1.*nm));
FILE *fppe = fopen(filenamer,"w+");
fprintf(fppe,"Total Angular momentum transfered (error), a: %.2gnm    v: %.2gc   Lmax: %d \n", a/(1.*nm), vv, b/(1.*nm), Lmax);
fprintf(fppe,"\n");
fprintf(fppe,"         betav       errDLEy                errDLHy                 errDLEsy                errDLHsy                 errDLEsyExt                errDLHsyExt\n");

for (int i = 0; i <= 9; ++i)
{
   //vv = (0.15 + (0.875-0.15)*i/17);
   vv = (0.5 + (0.95-0.5)*i/9);
   DL(r, vv, b, a, DLy);
   Ly=0.0;

   for (int rr = 0; rr < 6; ++rr){ 
      Ly += DLy[rr].real();}

   // Here print the total momentum
   fprintf(fpp,"%.17g %.17g %.17g %.17g %.17g  %.17g %.17g %.17g \n",vv,DLy[0].real(),DLy[1].real(),DLy[2].real(),DLy[3].real(),DLy[4].real(),DLy[5].real(), Ly);
   fprintf(fppe,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",vv,DLy[0].imag()/DLy[0].real(),DLy[1].imag()/DLy[1].real(),DLy[2].imag()/DLy[2].real(),DLy[3].imag()/DLy[3].real(),DLy[4].imag()/DLy[4].real(),DLy[5].imag()/DLy[5].real());
}}
// End the timer
auto end = std::chrono::high_resolution_clock::now();
// Calculate the duration
std::chrono::duration<double> duration = end - start;
// Print the duration in seconds
std::cout << "Execution time: " << duration.count() << " seconds." << std::endl;

cout << "a = " << a/nm << "nm." << endl;
cout << "r = " << r/nm << "nm." << endl;
cout << "v = " << vv << "c." << endl;
cout << "b = " << b/nm << "nm." << endl;
cout << endl;

/**for (int i = 0; i < 10; ++i)
{
    //r =(1.05 + i)*nm;
    //b = (1.5 + i)*nm;
    v = 0.1 + (0.99-0.1)*i/10;
    DL(r, v, b, a);
}**/

return(0);

}

//**************************************************** END ****************************************************
