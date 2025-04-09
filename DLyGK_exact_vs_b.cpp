//*****************************************************************************************************************
/*
DLGK_exact_vs_b.cpp

Description:
This program calculates the spectral angular momentum (dL/dω) transferred to a 
nanoparticle (NP) due to the interaction with a fast-moving electron. The calculation 
is performed by analytically solving the closed surface integral, resulting in a 
double multipolar sum expression. The program focuses on the external field 
contribution, the scattered field contribution and the interaction contribution.
This program calculates the total angular momentum transfer (DeltaL) for different 
values of the impact parameter b.

Key Features:
- The frequency integral is computed using the Gauss-Kronrod quadrature method, 
  with an adaptive partition for rapid convergence.
- The spectral results are saved to files named "dldw_*.dat".
- The integrated angular momentum results are saved to files named "DL_*.dat".

Usage:
- The user can input parameters manually or load them from a CSV file.
- The program supports multiple parameter sets and stops when the solution converges 
  to a specified relative error threshold or reaches the maximum multipolar order (Lmax = 50).
*/  
//	flags 
//      g++ -o OUT/DLyGK_exact_vs_b.out DLyGK_exact_vs_b.cpp -lcomplex_bessel -w 
//      g++ -o OUT/DLyGK_exact_vs_b.out DLyGK_exact_vs_b.cpp -lcomplex_bessel -w && ./OUT/DLyGK_exact_vs_b.out 
//
//          By Jesús Castrejon, jcastrejon@ciencias.unam.mx (25/02/2019)
// Modified by Jorge Luis Briseño, jorgeluisbrisenio@ciencias.unam.mx (21/02/2023)
//
//*****************************************************************************************************************
//*****************************************************************************************************************

#include "Irreducible_Integrals/IN51.h"
int Lmax = 51;

// Spherical Bessel (Hankel) funcions with *complex argument* (very important) for the scatterred fields. 
// Only found in SLATEC project (FORTRAN 77).

#include <complex_bessel_bits/sph_besselFunctions.h>  

// Source:       https://github.com/valandil/complex_bessel                                                                          

using namespace sp_bessel;
using namespace std;                  

typedef complex<double> dcomplex;          // Defines complex number with double entries.

double Cspeed = 137.035999139;
double Pi     = boost::math::constants::pi<double>();           // Parameter for Drude model of the NP and atomic units
double nm     = 100./5.2918;

//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************
#include "Dielectric_Functions/EpsilonDrudeAl.h"


//************************************************

static double IM[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IN[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IU[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IV[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IW[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IX[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IY[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IZ[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double ID[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double III[LSmax][LSmax+1][LSmax+1];


//**********************************************************************
// Common functions
//********************************************************************** 

//double BesselK(int n, double x){
//  return boost::math::cyl_bessel_k(n,x);
//}  

dcomplex BesselK(int n, double x){
    //return boost::math::cyl_bessel_k(n,x);
    dcomplex z = static_cast<dcomplex>(x);
    return besselK(n,z);
}                                                      

double LP(int l, int m, double x){
	return boost::math::legendre_p(l,m,x);
}

double factorial2(int n)
{
    if (n == 0 || n==1){
      return 1.;
    }else if(n > 1){
    return n*factorial2(n-2);
    }
    else{
        return pow(-1.,(n-1.)/2.)*n/factorial2(abs(n));
    }
}


double factorial(int n)
{
    if (n == 0 || n==1){
      return 1.;
    }else{
    return n*factorial(n-1);
    }
}


//********************************************************************************
// ****************** Gauss Konrod Quadrature weights and abscisas ********************
//********************************************************************************

//***********************************************************************
void gauss_konrod_w(double GKQt[2*nw1 + 1][3], double GKQf[2*nw2 + 1][3]){

  int const N1 = nw1;
  int const N2 = nw2;

   auto XGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::abscissa();    
   auto WGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::weights();      
   auto WGt =  boost::math::quadrature::gauss<double, N1>::weights();

   auto XGKf = boost::math::quadrature::gauss_kronrod<double, 2*N2+1>::abscissa();
   auto WGKf = boost::math::quadrature::gauss_kronrod<double, 2*N2+1>::weights();
   auto WGf =  boost::math::quadrature::gauss<double, N2>::weights();

   double WGLt[2*N1+1], WGLf[2*N2+1];


// Changes the n order array of Gaussian quadrature to a 2n+1 array

for(int i = 0; i < XGKt.size(); ++i){
      if (i % 2 == 0){
           WGLt[i] = WGt[i/2];}
      else WGLt[i] = 0.;
}

for(int i = 0; i < XGKf.size(); ++i){
      if (i % 2 == 0){
           WGLf[i] = WGf[i/2];}
      else WGLf[i] = 0.;
}

// writtes [0] abscisas, [1] Konrod weights, [2] Gauss weigths

for(int i = 0; i < 2*XGKt.size() - 1; ++i){
      if (i <= N1){ 
         GKQt[i][0] = -1.*XGKt[N1 - i];
         GKQt[i][1] = WGKt[N1 - i];                
         GKQt[i][2] = WGLt[N1 - i];
      }
      else{
         GKQt[i][0] = XGKt[-N1 + i];
         GKQt[i][1] = WGKt[-N1 + i];
         GKQt[i][2] = WGLt[-N1 + i];
      }
}

for(int i = 0; i < 2*XGKf.size() - 1; ++i){
      if (i <= N2){ 
         GKQf[i][0] = -1.*XGKf[N2 - i];
         GKQf[i][1] = WGKf[N2 - i];
         GKQf[i][2] = WGLf[N2 - i];
      }
      else{
         GKQf[i][0] = XGKf[-N2 + i];
         GKQf[i][1] = WGKf[-N2 + i];
         GKQf[i][2] = WGLf[-N2 + i];
      }
}

}





void gauss_konrod_w2(double GKQt[2*nw3 + 1][3], double GKQf[2*nw4 + 1][3]){

  int const N1 = nw3;
  int const N2 = nw4;

   auto XGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::abscissa();    
   auto WGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::weights();      
   auto WGt =  boost::math::quadrature::gauss<double, N1>::weights();

   auto XGKf = boost::math::quadrature::gauss_kronrod<double, 2*N2+1>::abscissa();
   auto WGKf = boost::math::quadrature::gauss_kronrod<double, 2*N2+1>::weights();
   auto WGf =  boost::math::quadrature::gauss<double, N2>::weights();

   double WGLt[2*N1+1], WGLf[2*N2+1];


// Changes the n order array of Gaussian quadrature to a 2n+1 array

for(int i = 0; i < XGKt.size(); ++i){
      if (i % 2 == 0){
           WGLt[i] = WGt[i/2];}
      else WGLt[i] = 0.;
}

for(int i = 0; i < XGKf.size(); ++i){
      if (i % 2 == 0){
           WGLf[i] = WGf[i/2];}
      else WGLf[i] = 0.;
}

// writtes [0] abscisas, [1] Konrod weights, [2] Gauss weigths

for(int i = 0; i < 2*XGKt.size() - 1; ++i){
      if (i <= N1){ 
         GKQt[i][0] = -1.*XGKt[N1 - i];
         GKQt[i][1] = WGKt[N1 - i];                
         GKQt[i][2] = WGLt[N1 - i];
      }
      else{
         GKQt[i][0] = XGKt[-N1 + i];
         GKQt[i][1] = WGKt[-N1 + i];
         GKQt[i][2] = WGLt[-N1 + i];
      }
}

for(int i = 0; i < 2*XGKf.size() - 1; ++i){
      if (i <= N2){ 
         GKQf[i][0] = -1.*XGKf[N2 - i];
         GKQf[i][1] = WGKf[N2 - i];
         GKQf[i][2] = WGLf[N2 - i];
      }
      else{
         GKQf[i][0] = XGKf[-N2 + i];
         GKQf[i][1] = WGKf[-N2 + i];
         GKQf[i][2] = WGLf[-N2 + i];
      }
}

}




//**************************************************************************************
// Recursive functions and parameter needed for calculate the scatterd fields
//*************************************************************************************

double II(int l, int m, int i1, int i2){

if (l == m - 2 || l == m - 1){return 0.;} 
     else if(l == m){ 	if (i2 % 2 == 0){
     		               return pow(-1.,m)*factorial2(2*m - 1)*boost::math::beta((i1 + m + 2)/2., (i2 + 1)/2.,1.);}
     		            else{return 0.;}}
     	else{return (1./(l - m))*((2.*l - 1.)*II(l - 1, m,i1, i2 + 1) - (l + m - 1.)*II(l - 2, m,i1, i2));}
              }


double alm(int l, int m){
	if (abs(m)<=l)
	{
		return pow((2.*l+1.)*factorial(l-m)/(4.*Pi*factorial(l+m)),0.5);
	}
	else{return 0.;}
}


dcomplex  A(int l, int m, double betav){

double gam = 1./sqrt(1.0 - pow(betav,2.0));
dcomplex  res = 0.;

if(m >= 0){
for (int j = m; j <= l; ++j){
	//res += pow(betav,-(l+1))*pow(1i,l-j)*factorial2(2*l+1)*alm(l,m)*III[l-1][m][j]
	///(pow(2.*gam,j)*factorial(l-j)*factorial((j-m)/2)*factorial((j+m)/2));
    res += pow(betav,-(l+1))*pow(1i,l-j)*factorial2(2*l+1)*III[l-1][m][j]
    /(pow(2.*gam,j)*factorial(l-j)*factorial((j-m)/2)*factorial((j+m)/2));
	}
return res;
 } 
else {
 	return pow(-1.,abs(m))*A(l,abs(m),betav);
 }
}



dcomplex  B(int l, int m, double betav){
return A(l,m+1,betav)*sqrt((l+m+1)*(l-m)) - A(l,m-1,betav)*sqrt((l-m+1)*(l+m));
}

dcomplex  PsiE(int l, int m, double w, double b,double betav, dcomplex BB){
double gam = 1./sqrt(1.0 - betav*betav);
return -2.*Pi*pow(1i,1-l)*w*BB*BesselK(m,w*b/(betav*Cspeed*gam))/(Cspeed*Cspeed*gam*l*(l+1));
}

dcomplex PsiM(int l, int m, double w, double b,double betav, dcomplex AA){
double gam = 1./sqrt(1.0 - betav*betav);
return -(4.*m)*Pi*pow(1i,1-l)*betav*w*AA*BesselK(m,w*b/(betav*Cspeed*gam))/(Cspeed*Cspeed*l*(l+1));
}


dcomplex j(int n, dcomplex z){    
return sph_besselJ(n,z);
}

dcomplex h(int n, dcomplex z){                              
return sph_hankelH1(n,z);
}


dcomplex djj(int n, dcomplex x){
	return j(n,x) + 0.5*x*(j(n-1,x)-j(n+1,x) - j(n,x)/x);
}

dcomplex dhh(int n, dcomplex x){
	return h(n,x) + 0.5*x*(h(n-1,x)-h(n+1,x) - h(n,x)/x);
}


// *********** Polarizabilities *********************
dcomplex tE(int l, double w, double a){
double x0 = w*a/Cspeed;
dcomplex xi = x0*sqrt(eps(w));
return -1i*( eps(w)*j(l, xi)*djj(l,x0) - j(l,x0)*djj(l,xi) )/( h(l,x0)*djj(l,xi) - eps(w)*j(l,xi)*dhh(l,x0) ); 
}

dcomplex tM(int l, double w, double a){
double x0 = w*a/Cspeed;
dcomplex xi = x0*sqrt(eps(w));
return -1i*( j(l, xi)*djj(l,x0) - j(l,x0)*djj(l,xi) )/( h(l,x0)*djj(l,xi) - j(l,xi)*dhh(l,x0) );
}



dcomplex fs(int l, dcomplex z){
double dl = 1.*l;
	return (dl + 1.)*1i*h(l,z)/z  - 1i*h(l+1,z);
}

dcomplex fe(int l, dcomplex z){
double dl = 1.*l;
  return (dl + 1.)*j(l,z)/z  - j(l+1,z);
}
//******************************************************************************************************
//******************************************************************************************************

void Omegas(double xi[2*NN + 4], double xk[2*NN + 4], double xg[2*NN + 4]){

double GKQw1[2*nw1 + 1][3], GKQw2[2*nw2 + 1][3];
gauss_konrod_w(GKQw1,GKQw2);
double GKQw3[2*nw3 + 1][3], GKQw4[2*nw4 + 1][3];
gauss_konrod_w2(GKQw3,GKQw4);

double w;

int NN1 = 2*nw1 + 1;
int NN2 = 2*nw2 + 1;
int NN3 = 2*nw3 + 1;
int NN4 = 2*nw4 + 1;

for(int i=0; i < NN1;  ++i){
xi[i] = ((w2-w1)/2.)*GKQw1[i][0] + ((w2+w1)/2.);;
xk[i] = ((w2-w1)/2.)*GKQw1[i][1];
xg[i] = ((w2-w1)/2.)*GKQw1[i][2];
}

for(int i=0; i < NN2;  ++i){
xi[i + NN1] = ((w3-w2)/2.)*GKQw2[i][0] + ((w3+w2)/2.);;
xk[i + NN1] = ((w3-w2)/2.)*GKQw2[i][1];
xg[i + NN1] = ((w3-w2)/2.)*GKQw2[i][2];
}

for(int i=0; i < NN3;  ++i){
xi[i + NN2+ NN1] = ((w4-w3)/2.)*GKQw3[i][0] + ((w4+w3)/2.);;
xk[i + NN2+ NN1] = ((w4-w3)/2.)*GKQw3[i][1];
xg[i + NN2+ NN1] = ((w4-w3)/2.)*GKQw3[i][2];
}

for(int i=0; i < NN4;  ++i){
xi[i + NN3+ NN1+ NN2] = ((w5-w4)/2.)*GKQw4[i][0] + ((w5+w4)/2.);;
xk[i + NN3+ NN1+ NN2] = ((w5-w4)/2.)*GKQw4[i][1];
xg[i + NN3+ NN1+ NN2] = ((w5-w4)/2.)*GKQw4[i][2];
}

}

void DL(double r, double vv, double b, double a, dcomplex DLy[6]){

double dm;
double dl;

double dldwy[6];
//dcomplex DLy[6];

dcomplex CM[Lmax][2*Lmax +1], tMl[Lmax];
dcomplex DE[Lmax][2*Lmax +1], tEl[Lmax];
dcomplex dz[2][Lmax], df[2][Lmax];

dcomplex IErty[4], IHrty[4];
dcomplex IErfy[4], IHrfy[4];

double dl1, dl2, dm1, dm2;
double normfactor; // A normalization factor due to the normalization of integrals

double IN1, IV1, IW1, IW2, IU1, IU2;
FUN(IM,IN,IU,IV,IW,IX,IY,IZ,ID,III);


// Calls the Gauss - Konrod function
double xi[2*NN + 4], xk[2*NN + 4], xg[2*NN + 4];
Omegas(xi, xk, xg);

double w, k0;

char filenamex[sizeof "Results/50/dldw_a1nm_v0.99c_b1.5nm_extF2.dat"];
sprintf(filenamex, "Results/%s/%d/dldw_a%.2gnm_v%.2g_b%.2gnm.dat", directoryLb, Lmax, a/(1.*nm), vv, b/(1.*nm));
FILE *fpx = fopen(filenamex,"w+");
fprintf(fpx,"Angular Momentum Spectrum, a: %.2gnm    v: %.2gc   b: %.2gnm    Lmax: %d \n", a/(1.*nm), vv, b/(1.*nm),Lmax);
fprintf(fpx,"\n");
fprintf(fpx,"         b         w(au)                   dlEdwy                   dlHdwy                 dlEsdwy                   dlHsdwy                   dlEsdwyExt                   dlHsdwyExt\n");

dcomplex AA[Lmax][2*Lmax+1], BB[Lmax][2*Lmax+1];
dcomplex aa, bb;

for (int l = 1; l <= Lmax; ++l){
    for (int m = -l; m <= l; ++m){ 
    AA[l-1][l+m] = A(l,m,vv);
    BB[l-1][l+m] = B(l,m,vv);
    }
}


for (int i = 0; i < 6; ++i){
DLy[i] = 0.;
}


for (int i = 0; i < iw4; ++i){     // Here we can select the Omega interval to integrate "2*NN + 4" for total.
    w = xi[i];
    k0 = w/Cspeed;

// Calculates  common parameters of the scatterred fields
for (int l = 1; l <= Lmax; l++){
    tEl[l-1] = tE(l,w,a);
    tMl[l-1] = tM(l,w,a);
    
    dz[0][l-1] = 1i*h(l,k0*r);
    df[0][l-1] = fs(l,k0*r);
    dz[1][l-1] = j(l,k0*r);
    df[1][l-1] = fe(l,k0*r);

    for (int m = -l; m <= l; ++m){

        aa = AA[l-1][m+l];
        bb = BB[l-1][l+m];
    
        DE[l-1][m + l] = pow(1i,l)*PsiE(l,m,w,b,vv, bb);
        CM[l-1][m + l] = pow(1i,l)*PsiM(l,m,w,b,vv, aa);

        //DE[l-1][m + l] = pow(1i,l)*alm(l,m)*PsiE(l,m,w,b,vv, bb);
        //CM[l-1][m + l] = pow(1i,l)*alm(l,m)*PsiM(l,m,w,b,vv, aa);
    
    }
} 



for (int rr = 0; rr < 6; ++rr){
dldwy[rr] = 0.;
}



for (int l1 = 1; l1 <= Lmax; l1++){
     dl1 = 1.*l1; 
for (int l2 = 1; l2 <= Lmax; l2++){
         dl2 = 1.*l2; 
         for (int m1 = -l1; m1 <= l1; m1++){
              dm1 = 1.*m1;
              for (int m2 = -l2; m2 <= l2; m2++){
                   dm2 = 1.*m2;

                    if(std::abs (m2-m1)==1){

                                    // Radial - Radial
                        normfactor = pow((2.*dl1+1.)*(l1+m1+1.)/((2.*dl1+3.)*(l1-m1+1.)),0.5);

                        IV1 = IV[l1-1][l2-1][m1+l1][m2+l2];

                        IW1 = IW[l1-1][l2-1][m1+l1][m2+l2];
                        //IW2 = IW[l1][l2-1][m1+l1+1][m2+l2]; 
                        IW2 = normfactor*IW[l1][l2-1][m1+l1+1][m2+l2];   

                        IU1 = IU[l1-1][l2-1][m1+l1][m2+l2];
                        //IU2 = IU[l1][l2-1][m1+l1+1][m2+l2];
                        IU2 = normfactor*IU[l1][l2-1][m1+l1+1][m2+l2];

                        /* Integrants, 0 -> Scat
                                       1 -> Ext
                                       2 -> Ext-Scat 
                                       3 -> Scat-Ext (second is conjugated)*/



                    // Zenital - Radial

                    // y component

                    IErty[0] = (- tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 
                              - (dl1-dm1+1.)*IU2) - tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);
     
                    IHrty[0] = (- tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1
                              - (dl1-dm1+1.)*IU2) + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);


                    IErty[1] =  0.0;
     
                    IHrty[1] =  0.0;


                    IErty[2] = (- tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 
                              - (dl1-dm1+1.)*IU2) - tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);
     
                    IHrty[2] = (- tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1
                              - (dl1-dm1+1.)*IU2) + tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);


                    IErty[3] = (- DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 
                              - (dl1-dm1+1.)*IU2) - CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);
     
                    IHrty[3] = (- CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1
                              - (dl1-dm1+1.)*IU2) + DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);


                    // Azimutal - Radial


                    // y compo

                    IErfy[0] =  -(dm1-dm2)*( tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 
                              - (dl1-dm1+1.)*IW2)+ tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1); 

                    IHrfy[0] =  -(dm1-dm2)*( -tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 
                              - (dl1-dm1+1.)*IW2)+ tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1); 

                    IErfy[1] = 0.0; 

                    IHrfy[1] =  0.0;  

                    IErfy[2] =  -(dm1-dm2)*( tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 
                              - (dl1-dm1+1.)*IW2)+ tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1); 

                    IHrfy[2] =  -(dm1-dm2)*( -tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 
                              - (dl1-dm1+1.)*IW2)+ tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);

                    IErfy[3] =  -(dm1-dm2)*( CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 
                              - (dl1-dm1+1.)*IW2)+ DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);  

                    IHrfy[3] =  -(dm1-dm2)*( -DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 
                              - (dl1-dm1+1.)*IW2)+ CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1); 

                    }


                else{
                  for (int rr = 0; rr < 4; ++rr){

                     IErty[rr] = 0.;
                     IHrty[rr] = 0.;

                     IErfy[rr] = 0.;
                     IHrfy[rr] = 0.;

                    }
                   }

                   /* Integrants (DL), 0 -> Scat
                                       1 -> Ext
                                       2 -> Ext-Scat 
                                       3 -> Scat-Ext (second is conjugated)*/


                        /* Integrants (dldw), 0 -> Ext-Scat (Electric)
                                       1 -> Ext-Scat (Magnetic)
                                       2 -> Scat-Scat (Electric) 
                                       3 -> Scat-Scat (Magnetic)
                                       4 -> Ext-Ext (Electric)
                                       5 -> Ext-Ext (Magnetic)*/

                    dldwy[0] +=  (1./(4.*Pi))*pow(r,3)*(( IErty[2] - IErfy[2] ).real()
                                                    +( IErty[3] - IErfy[3] ).real());
                    dldwy[1] +=  (1./(4.*Pi))*pow(r,3)*((IHrty[2] -IHrfy[2]).real()
                                                    +(IHrty[3] -IHrfy[3]).real());
                    dldwy[2] +=  (1./(4.*Pi))*pow(r,3)*( IErty[0] - IErfy[0] ).real();
                    dldwy[3] +=  (1./(4.*Pi))*pow(r,3)*(IHrty[0] -IHrfy[0]).real();
                    dldwy[4] +=  (1./(4.*Pi))*pow(r,3)*( IErty[1] - IErfy[1] ).real();
                    dldwy[5] +=  (1./(4.*Pi))*pow(r,3)*(IHrty[1] -IHrfy[1]).real();

            } // for m2
        }  // for m1

} // for l2
} //for l1

for (int rr = 0; rr < 6; ++rr){ 
  DLy[rr] += (xk[i]*(1. + 1i) - 1i*xg[i])*dldwy[rr];
}



// Here print the dldw's, for each, w, l
fprintf(fpx,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",b/(1.*nm),w,dldwy[0],dldwy[1],dldwy[2], dldwy[3], dldwy[4],dldwy[5]);

// cout << "In   " << i + 1 << "  of   " << 2*NN + 4 << endl;

}// for w

/* Integrants (dldw),                  0 -> Ext-Scat (Electric)
                                       1 -> Ext-Scat (Magnetic)
                                       2 -> Scat-Scat (Electric) 
                                       3 -> Scat-Scat (Magnetic)
                                       4 -> Ext-Ext (Electric)
                                       5 -> Ext-Ext (Magnetic)*/

//cout << "DLEy : " << DLy[0] << endl;
//cout << "DLHy : " << DLy[1] << endl;                // prints result
//cout << endl;

//cout << "DLEsy : " << DLy[2] << endl;
//cout << "DLHsy : " << DLy[3] << endl;                // prints result
//cout << endl;

cout << "DLy : " << DLy[0] + DLy[1] + DLy[2] + DLy[3] << endl;
//cout << endl;

// Close the files
fclose(fpx); 

} //end void



//*********************************************************************************************************
//********************************************** MAIN PROGRAM *********************************************
int main(void){

cout.precision(17);

int choice;

cout << endl;
cout << R"(DLGK_exact_vs_b.cpp

Description:
This program calculates the spectral angular momentum (dL/dω) transferred to a 
nanoparticle (NP) due to the interaction with a fast-moving electron. The calculation 
is performed by analytically solving the closed surface integral, resulting in a 
double multipolar sum expression. The program focuses exclusively on the external field 
contribution.This program calculates the total angular momentum transfer (DeltaL) for different 
values of the impact parameter b.

Key Features:
- The frequency integral is computed using the Gauss-Kronrod quadrature method, 
  with an adaptive partition for rapid convergence.
- The spectral results are saved to files named "dldw_*.dat".
- The integrated angular momentum results are saved to files named "DL_*.dat".

Usage:
- The user can input parameters manually or load them from a CSV file.
- The program supports multiple parameter sets and stops when the solution converges 
  to a specified relative error threshold or reaches the maximum multipolar order (Lmax = 50).

  )" << endl << endl;

cout << "Select an option to input parameters:" << endl;
    cout << R"(1. Select the following parameters
                double a = 5.0*nm;
                double bInit = 5.5*nm;
                double bFin = 10.0*nm;
                double r = 5.05*nm;
                double vv = 0.7;)" << endl;
    cout << R"(2. Select the following parameters
                double a = 10.0*nm; 
                double bInit = 11.5*nm;
                double bFin = 20.0*nm;
                double r = 10.05*nm;
                double vv = 0.7;)" << endl;
    cout << R"(3. Select the following parameters
                double a = 20.0*nm;
                double bInit = 21.5*nm;
                double bFin = 30.0*nm;                             
                double r = 20.05*nm;
                double vv = 0.7;)" << endl;
    cout << R"(4. Select the following parameters
                double a = 50.0*nm;
                double bInit = 51.5*nm;
                double bFin = 60.0*nm;                             
                double r = 50.05*nm;
                double vv = 0.7;)" << endl;
    cout << R"(5. Load parameters from parameters_vs_b.csv file
                Note that the file should be written in the form:
                    a,bInit,bFin,r,vv
                    50.0,51.5,60.0,50.05,0.5)" << endl;
    cout << "Option:  ";
    cin >> choice;
cout << endl;

double a, b, bInit, bFin, r, vv;

if (choice == 1) {                        
        a = 5.0*nm;     
        bInit = 5.5*nm;
        bFin = 10.0*nm;                   
        r = 5.05*nm;
        vv = 0.7;
    } 
else if (choice == 2) {                        
        a = 10.0*nm;
        bInit = 11.5*nm;
        bFin = 20.0*nm;                        
        r = 10.05*nm;
        vv = 0.7;
}
else if (choice == 3) {                        
        a = 20.0*nm;
        bInit = 21.5*nm;
        bFin = 30.0*nm;                        
        r = 20.05*nm;
        vv = 0.7;
}
else if (choice == 4) {                        
        a = 50.0*nm;
        bInit = 51.5*nm;
        bFin = 60.0*nm;                       
        r = 50.05*nm;
        vv = 0.7;
}
else if (choice == 5) {
        string filename = "Parameters/parameters_vs_b.csv";
        cout << R"(Reading the parameters in the parameters_vs_b.csv filename. 
                    Note that the file should be written in the form:
                    a,bInit, bFin,r,vv
                    50.0,51.5,60.0,50.05,0.5)";
        //cin >> filename;

        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file!" << endl;
            return 1;
        }

        string line;
        getline(file, line); // Skip header
        getline(file, line);

        stringstream ss(line);
        string value;

        getline(ss, value, ',');
        cout << endl;
        cout << endl;
        cout << "Parameters: \n \n";
        cout << "a = " << stod(value) << " nm\n";
        a = stod(value) * nm;

        getline(ss, value, ',');
        cout << "bInit = " << stod(value) << " nm \n";
        bInit = stod(value) * nm;

        getline(ss, value, ',');
        cout << "bFin = " << stod(value) << " nm \n";
        bFin = stod(value) * nm;

        getline(ss, value, ',');
        cout << "r = " << stod(value) << " nm \n";
        r = stod(value) * nm;

        getline(ss, value, ',');
        cout << "vv = " << stod(value) << "*c \n";
        vv = stod(value);

        cout << endl;
        cout << endl;
        file.close();
    } else {
        cerr << "Invalid choice!" << endl;
        return 1;
    }

dcomplex DLy[6];
const double errorThreshold = 0.0001; // 0.01% threshold
const int maxLmax = 51; // Maximum multipolar order
const int numOfPoints = 9; // Number of points for b
int minLmax; // Minimum multipolar order

cout << "Select minimum value of multipolar order L (suggested-->1):  ";
    cin >> minLmax;
cout << endl;

// Vector to store Ly values for each i in the previous Lmax iteration
vector<double> previousLy(numOfPoints + 1, 0.0); // Initialize to 0.0

char filenamel[sizeof "Results/MultipolarConvergence_a1nm_v0.99c_b1.5nm_error_extF2.dat"];
sprintf(filenamel, "Results/%s/MultipolarConvergence_a%.2gnm_v%.2g.dat", directoryLb, a/(1.*nm), vv);
FILE *fppl = fopen(filenamel,"w+");
fprintf(fppl,"Multipolar Convergence analysis, a: %.2gnm    v: %.2gc", a/(1.*nm), vv);
fprintf(fppl,"\n");
fprintf(fppl,"         Lmax       Average Relative Error\n");

// Start the timer
auto start = std::chrono::high_resolution_clock::now();
//
// Start loop for Lmax
//
for (Lmax = minLmax; Lmax <= maxLmax; Lmax++){
//cout << "Lmax = " << Lmax << endl;
// Print dynamic line
cout << "Progress: " 
          << "Lmax : " << Lmax
          << " | v : " << vv << "*c"
          << ", NP radius: " << a/nm << " nm" << endl;
char filename[sizeof "Results/50/DL_a1nm_v0.99c_b1.5nm_extF2.dat"];
sprintf(filename, "Results/%s/%d/DL_a%.2gnm_v%.2g.dat", directoryLb, Lmax,a/(1.*nm), vv);
FILE *fpp = fopen(filename,"w+");
fprintf(fpp,"Total Angular momentum transfered, a: %.2gnm    v: %.2gc   Lmax: %d  \n", a/(1.*nm), vv, Lmax);
fprintf(fpp,"\n");
fprintf(fpp,"         b         DLEy                  DLHy                  DLEsy                  DLHsy                  DLEsyExt                  DLHsyExt                  DLy\n");


char filenamer[sizeof "Results/50/DL_a1nm_v0.99c_b1.5nm_error_extF2.dat"];
sprintf(filenamer, "Results/%s/%d/DL_a%.2gnm_v%.2g_error.dat", directoryLb, Lmax, a/(1.*nm), vv);
FILE *fppe = fopen(filenamer,"w+");
fprintf(fppe,"Total Angular momentum transfered (error), a: %.2gnm    v: %.2gc   Lmax: %d \n", a/(1.*nm), vv, Lmax);
fprintf(fppe,"\n");
fprintf(fppe,"         b       errDLEy                errDLHy                 errDLEsy                errDLHsy                 errDLEsyExt                errDLHsyExt\n");

double totalRelativeError = 0.0; // Accumulator for relative errors
bool thresholdReached = false;

for (int i = 0; i <= numOfPoints; ++i)
{   
    b = (bInit + (bFin-bInit)*i/numOfPoints);

   DL(r, vv, b, a, DLy);
   double Ly=0.0;

   for (int rr = 0; rr < 6; ++rr){ 
   Ly += DLy[rr].real();}

   // Calculate the error
    double errorLy = fabs(Ly - previousLy[i]);
    double relativeError = errorLy / fabs(Ly);

    // Accumulate the relative error
    totalRelativeError += relativeError;

   // Here print the total momentum
   fprintf(fpp,"%.17g %.17g %.17g %.17g %.17g  %.17g %.17g %.17g \n",b/(1.*nm),DLy[0].real(),DLy[1].real(),DLy[2].real(),DLy[3].real(),DLy[4].real(),DLy[5].real(), Ly);
   fprintf(fppe,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",b/(1.*nm),DLy[0].imag(),DLy[1].imag(),DLy[2].imag(),DLy[3].imag(),DLy[4].imag(),DLy[5].imag());

   // Update previousLy for the next iteration
    previousLy[i] = Ly;
        }

        // Calculate the average relative error
        double averageRelativeError = totalRelativeError / (numOfPoints + 1);
        fprintf(fppl,"%d %.17g\n",Lmax,averageRelativeError);
        cout << "Average Relative Error is: " << averageRelativeError << endl;
        cout << "------------------------------------------------------" << endl << endl;

        // Check if the average error is below the threshold
        if (averageRelativeError < errorThreshold) {
            thresholdReached = true;
            fprintf(fppl,"Convergence achieved at Lmax =%d with average relative error = %.17g\n",Lmax,averageRelativeError);
            cout << "Convergence achieved at Lmax = " << Lmax << " with average relative error = " << averageRelativeError << endl;
            cout << "------------------------------------------------------" << endl << endl;
        }

        // Close the files
        fclose(fpp);
        fclose(fppe);

        // Check if the threshold was reached
        if (thresholdReached) {
            break; // Exit the outer loop
        }
    }


// End the timer
auto end = std::chrono::high_resolution_clock::now();
// Calculate the duration
std::chrono::duration<double> duration = end - start;
// Print the duration in seconds
std::cout << "Execution time: " << duration.count() << " seconds." << std::endl;
fprintf(fppl,"Execution time: %.17g seconds.\n",duration.count());
fclose(fppl);

cout << "a = " << a / nm << "nm." << endl;
cout << "r = " << r / nm << "nm." << endl;
cout << "v = " << vv << "c." << endl;
cout << "bInit = " << bInit / nm << "nm." << endl;
cout << "bFin = " << bFin / nm << "nm." << endl;
cout << endl;
 
return(0);

}

//**************************************************** END ****************************************************
