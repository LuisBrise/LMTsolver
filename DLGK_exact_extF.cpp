//*****************************************************************************************************************
//
//   DPGK_exact.cpp:
//
//    This program calculates the spectral angular momentum dL/dw transferred to a 
//   NP due a fast swfit electron. It calculate the closed surface integral 
//   analitically, obteining a double multipolar sum espresion. The spectrum is 
//   written to a file named "dLdw*.dat". The program ONLY considers the external field.
//   Also, the integral in frequencies is calculated via Gauss - Kronrod using
//   a partition set to fast convergence. The results are written to a file named "DL*.dat". 
//  
//	flags 
//      g++ -o DLGK_exact_extF.out DLGK_exact_extF.cpp -lcomplex_bessel -w
//
//          By Jesús Castrejon, jcastrejon@ciencias.unam.mx (25/02/2019)
// Modified by Jorge Luis Briseño, jorgeluisbrisenio@ciencias.unam.mx (22/02/2023)
//
//*****************************************************************************************************************
//*****************************************************************************************************************


#include "IN11.h"

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
double wp     = 0.555;
double Gamma  = 0.00555;



//************************************************
// Maximum number of multipole to take into account

const int Lmax = 5;

// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|
const int nw1 = 51;   
const int nw2 = 201;   
const int nw3 = 51;                   // Drude Al
const int nw4 = 51;  
const int NN = nw1 + nw2 + nw3 + nw4;

double w1 = 0.;
double w2 = .3;
double w3 = .4;
double w4 = 2.;
double w5 = 70.; // Omega cut for freq. integral

int iw1 = 2*nw1 + 1;
int iw2 = 2*(nw1 + nw2) + 2;
int iw3 = 2*(nw1 + nw2 +nw3) + 3;
int iw4 = 2*NN + 4;


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
// *********************** Dielectric function ***************************
//**********************************************************************

dcomplex eps(double w){
return 1. - pow(wp,2.)/(w*(w + 1i*Gamma));
} 


//**********************************************************************
// Common functions
//********************************************************************** 

double BesselK(int n, double x){
	return boost::math::cyl_bessel_k(n,x);
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
	res += pow(betav,-(l+1))*pow(1i,l-j)*factorial2(2*l+1)*alm(l,m)*III[l-1][m][j]
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




void DL(double r, double vv, double b, double a){

double dm;
double dl;

double dLdwx[6], dLdwy[6], dLdwz[6];
dcomplex DLx[6], DLy[6], DLz[6];

dcomplex CM[Lmax][2*Lmax +1], tMl[Lmax];
dcomplex DE[Lmax][2*Lmax +1], tEl[Lmax];
dcomplex dz[2][Lmax], df[2][Lmax];

dcomplex IEStr[4], IESfr[4], IHStr[4], IHSfr[4];
dcomplex IECtr[4], IHCtr[4];
dcomplex IECSfr[4], IHCSfr[4];
dcomplex IECCfr[4], IHCCfr[4];

double dl1, dl2, dm1, dm2;

double IN1, IV1, IW1, IW2, IW3, IU1, IU2, IU3, IU4;
double IM1, IZ1, IX1, IX2, IX3, IY1, IY2, IY3, IY4, ID1, ID2;
FUN(IM,IN,IU,IV,IW,IX,IY,IZ,ID,III);


// Calls the Gauss - Konrod function
double xi[2*NN + 4], xk[2*NN + 4], xg[2*NN + 4];
Omegas(xi, xk, xg);

double w, k0;


char filename[sizeof "Results/DL_a1nm_v0.99c_b1.5nm_extF.dat"];
sprintf(filename, "DL_a%.2gnm_v%.2g_b%.2gnm_extF.dat", a/(1.*nm), vv, b/(1.*nm));
FILE *fpp = fopen(filename,"w+");
fprintf(fpp,"Total angular momentum transfered, a: %.2gnm    v: %.2gc   b: %.2gnm   Lmax: %d  \n", a/(1.*nm), vv, b/(1.*nm), Lmax);
fprintf(fpp,"\n");
fprintf(fpp,"         DLEx                    DLHx                  DLEz                  DLHz            DLEsx                   DLHsx                 DLEsz                 DLHsz\n");


char filenamer[sizeof "Results/DL_a1nm_v0.99c_b1.5nm_error_extF.dat"];
sprintf(filenamer, "DL_a%.2gnm_v%.2g_b%.2gnm_error_extF.dat", a/(1.*nm), vv, b/(1.*nm));
FILE *fppe = fopen(filenamer,"w+");
fprintf(fppe,"Total angular momentum transfered, a: %.2gnm    v: %.2gc   b: %.2gnm      Lmax: %d \n", a/(1.*nm), vv, b/(1.*nm), Lmax);
fprintf(fppe,"\n");
fprintf(fppe,"       errDLEx                   errDLHx                 errDLEz                errDLHz        errDLEsx               errDLHsx               errDLEsz              errDLHsz\n");


char filenamex[sizeof "Results/dLdw_a1nm_v0.99c_b1.5nm_extF.dat"];
sprintf(filenamex, "dLdw_a%.2gnm_v%.2g_b%.2gnm_extF.dat", a/(1.*nm), vv, b/(1.*nm));
FILE *fpx = fopen(filenamex,"w+");
fprintf(fpx,"Angular Momentum Spectrum, a: %.2gnm    v: %.2gc   b: %.2gnm    Lmax: %d \n", a/(1.*nm), vv, b/(1.*nm),Lmax);
fprintf(fpx,"\n");
fprintf(fpx,"         w(au)                   dLEdwx                 dLHdwx                 dLEdwz                dLHdwz                  dLEsdwx                dLHsdwx                dLEsdwz                dLHsdwz              dLEedwx                dLHedwx                dLEedwz                dLHedwz\n");



dcomplex AA[Lmax][2*Lmax+1], BB[Lmax][2*Lmax+1];
dcomplex aa, bb;

for (int l = 1; l <= Lmax; ++l){
    for (int m = -l; m <= l; ++m){ 
    AA[l-1][l+m] = A(l,m,vv);
    BB[l-1][l+m] = B(l,m,vv);
    }
}


for (int i = 0; i < 6; ++i){
DLx[i] = 0.;
DLy[i] = 0.;
DLz[i] = 0.;
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
    
        DE[l-1][m + l] = pow(1i,l)*alm(l,m)*PsiE(l,m,w,b,vv, bb);
        CM[l-1][m + l] = pow(1i,l)*alm(l,m)*PsiM(l,m,w,b,vv, aa);
    
    }
} 



for (int rr = 0; rr < 6; ++rr){
dLdwx[rr] = 0.;
dLdwy[rr] = 0.;
dLdwz[rr] = 0.;
}

for (int l1 = 1; l1 <= Lmax; ++l1){
     dl1 = 1.*l1; 
for (int l2 = 1; l2 <= Lmax; ++l2){
         dl2 = 1.*l2; 
         for (int m1 = -l1; m1 <= l1; ++m1){
              dm1 = 1.*m1;
              for (int m2 = -l2; m2 <= l2; ++m2){
                   dm2 = 1.*m2;


                    if(m2 == m1+1 || m2 == m1-1){

                                    // Radial - Radial 

                        //IN1 = IN[l1-1][l2-1][m1+l1][m2+l2];
                        IV1 = IV[l1-1][l2-1][m1+l1][m2+l2];

                        IW1 = IW[l1-1][l2-1][m1+l1][m2+l2];
                        IW2 = IW[l1][l2-1][m1+l1+1][m2+l2];   
                        //IW3 = IW[l1-1][l2][m1+l1][m2+l2+1];

                        IU1 = IU[l1-1][l2-1][m1+l1][m2+l2];
                        IU2 = IU[l1][l2-1][m1+l1+1][m2+l2];
                        //IU3 = IU[l1-1][l2][m1+l1][m2+l2+1];
                        //IU4 = IU[l1][l2][m1+l1+1][m2+l2+1];

                        /* Integrants, 0 -> Scat
                                       1 -> Ext
                                       2 -> Ext-Scat 
                                       3 -> Scat-Ext (second is conjugated)*/


                    IEStr[0] = 0; 

                    IHStr[0] = 0; 


                    IEStr[1] = 1i*(dm1-dm2)*r*r*dl2*(dl2+1.)*( -CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0))*dm1*IU1 
                                - DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0))*( (dl1+1.)*IW1 -(dl1-dm1+1.)*IU2 ) ); 

                    IHStr[1] = 1i*(dm1-dm2)*r*r*dl2*(dl2+1.)*( +DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0))*dm1*IU1 
                                - CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0))*( (dl1+1.)*IW1 -(dl1-dm1+1.)*IU2 ) );



                    IEStr[2] = 0; //G

                    IHStr[2] = 0; //G



                    IEStr[3] = 0; //G

                    IHStr[3] = 0; //G


                                   // Zenital - Radial

                    IECtr[0] = 0;
     
                    IHCtr[0] =  0;



                    IECtr[1] = r*r*dl2*(dl2+1.)*( -CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0))*dm1*IU1 
                                - DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0))*( (dl1+1.)*IW1 -(dl1-dm1+1.)*IU2 ) ); 
     
                    IHCtr[1] =  r*r*dl2*(dl2+1.)*( +DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0))*dm1*IU1 
                                - CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0))*( (dl1+1.)*IW1 -(dl1-dm1+1.)*IU2 ) );



                    IECtr[2] = 0;
     
                    IHCtr[2] = 0;




                    IECtr[3] = 0;
     
                    IHCtr[3] = 0;




                                    // Azimutal - Radial


                    IECSfr[0] =  0; 

                    IHCSfr[0] =  0;  

                    IECSfr[1] =  -(dm1-dm2)*r*r*dl2*(dl2+1.)*( +CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0))*( (dl1+1.)*IV1 - (dl1-dm1+1.)*IW2 ) 
                                + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0))*IW1 );

                    IHCSfr[1] =  -(dm1-dm2)*r*r*dl2*(dl2+1.)*( -DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0))*( (dl1+1.)*IV1 - (dl1-dm1+1.)*IW2 ) 
                                + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0))*IW1 ); 

                    IECSfr[2] =  0; 

                    IHCSfr[2] =  0;

                    IECSfr[3] =  0; 

                    IHCSfr[3] =  0;





                                        //  Zenital - Radial

                    IECCfr[0] =  0; 

                    IHCCfr[0] =  0; 



                    IECCfr[1] =  1i*r*r*dl2*(dl2+1.)*( +CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0))*( (dl1+1.)*IV1 - (dl1-dm1+1.)*IW2 ) 
                                + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0))*IW1 );

                    IHCCfr[1] =  1i*r*r*dl2*(dl2+1.)*( -DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0))*( (dl1+1.)*IV1 - (dl1-dm1+1.)*IW2 ) 
                                + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0))*IW1 ); 


                    IECCfr[2] =  0; 

                    IHCCfr[2] =  0; 


                    IECCfr[3] =  0; 

                    IHCCfr[3] =  0; }


                else{
                  for (int rr = 0; rr < 4; ++rr){
                     IEStr[rr] = 0.;
                     IHStr[rr] = 0.;
                     IECtr[rr] = 0.;
                     IHCtr[rr] = 0.;
                     IECSfr[rr] = 0.;
                     IHCSfr[rr] = 0.;
                     IECCfr[rr] = 0.;
                     IHCCfr[rr] = 0.;
                    }
                   }




                if(m2 == m1){
                                              // Radial - Radial 

                        //IM1 = IM[l1-1][l2-1][m1+l1][m2+l2];
                        //IZ1 = IZ[l1-1][l2-1][m1+l1][m2+l2];

                        ID1 = ID[l1-1][l2-1][m1+l1][m2+l2];
                        ID2 = ID[l1][l2-1][m1+l1+1][m2+l2];

                        //IX1 = IX[l1-1][l2-1][m1+l1][m2+l2];
                        //IX2 = IX[l1][l2-1][m1+l1+1][m2+l2];
                        //IX3 = IX[l1-1][l2][m1+l1][m2+l2+1];

                        //IY1 = IY[l1-1][l2-1][m1+l1][m2+l2];
                        //IY2 = IY[l1][l2-1][m1+l1+1][m2+l2];
                        //IY3 = IY[l1-1][l2][m1+l1][m2+l2+1];
                        //IY4 = IY[l1][l2][m1+l1+1][m2+l2+1];



                    IESfr[0] = 0; // G of "Good" cheked with TestGKL.cpp

                    IHSfr[0] = 0; //G


                    IESfr[1] = 2.*1i*r*r*dl2*(dl2+1.)*( +CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0))*( (dl1+1.)*IM1 - (dl1-dm1+1.)*ID2 ) 
                                + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0))*ID1 ); // G of "Good" cheked with TestGKL.cpp

                    IHSfr[1] = 2.*1i*r*r*dl2*(dl2+1.)*( -DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0))*( (dl1+1.)*IM1 - (dl1-dm1+1.)*ID2 ) 
                                + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0))*ID1 ); //G


                    IESfr[2] = 0; // G of "Good" cheked with TestGKL.cpp

                    IHSfr[2] = 0; //G


                    IESfr[3] = 0; // G of "Good" cheked with TestGKL.cpp

                    IHSfr[3] = 0; //G 
                }  
                                              
                else{
                  for (int rr = 0; rr < 4; ++rr){
                     IESfr[rr] = 0.;
                     IHSfr[rr] = 0.;
                    }
                   }




                     
                    dLdwx[0] += (1./(4.*Pi))*(1i).real();
                    dLdwz[0] += (1./(4.*Pi))*(1i).real();

                    dLdwx[1] +=  (1./(4.*Pi))*(1i).real();
                    dLdwz[1] +=  (1./(4.*Pi))*(1i).real();

                    dLdwx[2] +=  (1./(4.*Pi))*(1i).real();
                    dLdwz[2] +=  (1./(4.*Pi))*(1i).real();

                    dLdwx[3] +=  (1./(4.*Pi))*(1i).real();
                    dLdwz[3] +=  (1./(4.*Pi))*(1i).real();

                    // ONLY THE FOLLOWING DLDW ARE CORRECT --------------->

                    dLdwx[4] +=  (1./(4.*Pi))*( IEStr[1] + IECCfr[1] ).real();
                    dLdwy[4] +=  (1./(4.*Pi))*( IECtr[1] - IECSfr[1] ).real();
                    dLdwz[4] +=  (1./(4.*Pi))*IESfr[1].real();

                    dLdwx[5] +=  (1./(4.*Pi))*( IHStr[1] + IHCCfr[1] ).real();
                    dLdwy[5] +=  (1./(4.*Pi))*( IHCtr[1] - IHCSfr[1] ).real();
                    dLdwz[5] +=  (1./(4.*Pi))*( IHSfr[1] ).real();

                    //cout << "DLx = " << dLdwx[4] << endl;


            } // for m2
        }  // for m1

} // for l2
} //for l1




for (int rr = 0; rr < 6; ++rr){ 
  DLx[rr] += (xk[i]*(1. + 1i) - 1i*xg[i])*dLdwx[rr];
  DLy[rr] += (xk[i]*(1. + 1i) - 1i*xg[i])*dLdwy[rr];
  DLz[rr] += (xk[i]*(1. + 1i) - 1i*xg[i])*dLdwz[rr];
}




// Here print the dpdw's, for each, w, l
fprintf(fpx,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",w,dLdwx[0],dLdwx[1],dLdwz[0],dLdwz[1], dLdwx[2], dLdwx[3], dLdwz[2], dLdwz[3],dLdwx[4],dLdwx[5],dLdwz[4],dLdwz[5]);

//cout << "In   " << i + 1 << "  of   " << 2*NN + 4 << endl;
//cout << "I ext = " << IHCCfr[1].real() << "  " << IECCfr[1].real() << endl;


}// for w

// Here print the total momentum
fprintf(fpp,"%.17g %.17g %.17g %.17g  %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",DLx[0].real(),DLx[1].real(),DLz[0].real(),DLz[1].real(), DLx[2].real(),DLx[3].real(),DLz[2].real(),DLz[3].real(),DLx[4].real(),DLx[5].real(),DLz[4].real(),DLz[5].real());
fprintf(fppe,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",DLx[0].imag(),DLx[1].imag(),DLz[0].imag(),DLz[1].imag(), DLx[2].imag(),DLx[3].imag(),DLz[2].imag(),DLz[3].imag(),DLx[4].imag(),DLx[5].imag(),DLz[4].imag(),DLz[5].imag());


cout << endl;
//cout << "DLEx : " << DLx[0] << endl;
//cout << "DLHx : " << DLx[1] << endl;                // prints result
//cout << endl;
//cout << "DLEz : " << DLz[0] << endl;
//cout << "DLHz : " << DLz[1] << endl;
//cout << endl;
//cout << "DLEsx : " << DLx[2] << endl;
//cout << "DLHsx : " << DLx[3] << endl;                // prints result
//cout << endl;
//cout << "DLEsz : " << DLz[2] << endl;
//cout << "DLHsz : " << DLz[3] << endl;
//cout << endl;
cout << "DLEex : " << DLx[4] << endl;
cout << "DLHex : " << DLx[5] << endl;                // prints result
cout << endl;
cout << "DLEey : " << DLy[4] << endl;
cout << "DLHey : " << DLy[5] << endl;                // prints result
cout << endl;
cout << "DLEez : " << DLz[4] << endl;
cout << "DLHez : " << DLz[5] << endl;
cout << endl;
cout << endl;

} //end void



//*********************************************************************************************************
//********************************************** MAIN PROGRAM *********************************************
int main(void){


cout.precision(17);

cout << endl;
cout << "Momentum Gauss - Kronrod :" << endl;
cout << endl;
cout << "Lmax = " << Lmax << endl;
cout << endl;

double b = 2.5*nm;                        
double a = 1.*nm;                        
double r = 1.05*nm;
double v = 0.5;




cout << "a = " << a/nm << "nm." << endl;
cout << "v = " << v << "c." << endl;
cout << "b = " << b/nm << "nm." << endl;
cout << endl;

DL(r, v, b, a); 


 
return(0);


}

//**************************************************** END ****************************************************
