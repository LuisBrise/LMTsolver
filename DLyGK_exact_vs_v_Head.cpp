//*****************************************************************************************************************
//
/*
DLGK_exact_vs_v_Head.cpp

Description:
This program calculates the spectral angular momentum (dL/dω) transferred to a 
nanoparticle (NP) due to the interaction with a fast-moving electron. The calculation 
is performed by analytically solving the closed surface integral, resulting in a 
double multipolar sum expression. The program focuses on the external field 
contribution, the scattered field contribution and the interaction contribution.
This program calculates the total angular momentum transfer (DeltaL) for different 
values of the velocity of the electron v.

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
//  flags 
//      g++ -o OUT/quad_DLyGK_exact_vs_v.out DLyGK_exact_vs_v_Head.cpp -lcomplex_bessel -w 
//      g++ -o OUT/DLyGK_exact_vs_v.out DLyGK_exact_vs_v_Head.cpp -lcomplex_bessel -w && ./OUT/DLyGK_exact_vs_v.out 
//
//          By Jesús Castrejon, jcastrejon@ciencias.unam.mx (25/02/2019)
// Modified by Jorge Luis Briseño, jorgeluisbrisenio@ciencias.unam.mx (21/02/2023)
//
//*****************************************************************************************************************
//*****************************************************************************************************************
#include "Irreducible_Integrals/IN51.h"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

using namespace std;  


    using std::exp;
    using std::log;
    using std::pow;
    using std::sqrt;
    using std::sin;
    using std::cos;
    using std::conj;


typedef boost::multiprecision::cpp_complex<50> dcomplex;          // Defines complex number with double entries.
typedef boost::multiprecision::cpp_bin_float_50  dfloat;


//typedef boost::multiprecision::cpp_complex_quad dcomplex;
//typedef boost::multiprecision::cpp_bin_float_quad dfloat;


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

//const int Lmax = 3;

//************************************************
dfloat Cspeed = 137.035999084;
dfloat Pi     = boost::math::constants::pi<dfloat>();           // Parameter for Drude model of the NP and atomic units
dfloat nm     = dfloat(100)/5.2918;
//dfloat wp     = 0.555;
//dfloat Gamma  = 0.00555;
dcomplex imag_unit{0, 1};


//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************
#include "Dielectric_Functions_dfloat/EpsilonDrudeAl.h"


dfloat BesselK(int n, dfloat x){
    return boost::math::cyl_bessel_k(n,x);
}                                                  



dfloat LP(int l, int m, dfloat x){
  return boost::math::legendre_p(l,m,x);
}

dfloat factorial2(int n){
   return boost::math::double_factorial<dfloat>(n);
}

dfloat factorial(int n){
  return  boost::math::factorial<dfloat>(n);
}



//**************************************************************************************
// Recursive functions and parameter needed for calculate the scatterd fields
//*************************************************************************************

dfloat alm(int l, int m){
  if (abs(m)<=l)
  {
    return pow((2.*l+1.)*factorial(l-m)/(4.*Pi*factorial(l+m)),0.5);
  }
  else{return 0.;}
}


dcomplex  A(int l, int m, dfloat betav){

dfloat gam = 1./sqrt(1.0 - pow(betav,2.0));
dcomplex  res = 0.;

if(m >= 0){
for (int j = m; j <= l; ++j){
  res += pow(betav,-(l+1))*pow(imag_unit,l-j)*factorial2(2*l+1)*alm(l,m)*III[l-1][m][j]
  /(pow(2.*gam,j)*factorial(l-j)*factorial((j-m)/2)*factorial((j+m)/2));
  }
return res;
 } 
else {
  return pow(-1.,abs(m))*A(l,abs(m),betav);
 }
}

dcomplex  B(int l, int m, dfloat betav){
return A(l,m+1,betav)*sqrt((l+m+1)*(l-m)) - A(l,m-1,betav)*sqrt((l-m+1)*(l+m));
}

dcomplex  PsiE(int l, int m, dfloat w, dfloat b,dfloat betav, dcomplex BB){
dfloat gam = dfloat(1)/sqrt(dfloat(1) - betav*betav);
return -dfloat(2)*Pi*pow(imag_unit,1-l)*w*BB*BesselK(m,w*b/(betav*Cspeed*gam))/(Cspeed*Cspeed*gam*l*(l+1));
}

dcomplex PsiM(int l, int m, dfloat w, dfloat b,dfloat betav, dcomplex AA){
dfloat gam = dfloat(1)/sqrt(dfloat(1) - betav*betav);
return -(dfloat(4)*m)*Pi*pow(imag_unit,1-l)*betav*w*AA*BesselK(m,w*b/(betav*Cspeed*gam))/(Cspeed*Cspeed*l*(l+1));
}



complex<double> double_j(int n, complex<double>  z){   
if( z == 0. ){
    if(n==0){ return 1.;}
      else{ return 0.; }
}
 else{return sp_bessel::sph_besselJ(n,z);}
}


complex<double> double_h(int n, complex<double>  z){   
if( z == 0. ){
      if(n==0){ return 1.;}
          else{ return 0.; }}
 else{return sp_bessel::sph_hankelH1(n,z);}
}


dcomplex j(int n, dcomplex  zz){  
   complex<double> z = static_cast<complex<double>>(zz);
   complex<double> jj =  double_j(n, z);
   return static_cast<dcomplex>(jj);
}


dcomplex h(int n, dcomplex  zz){  
   complex<double> z = static_cast<complex<double>>(zz);
   complex<double> hh =  double_h(n, z);
   return static_cast<dcomplex>(hh);
}


dcomplex djj(int n, dcomplex x){
  return j(n,x) + 0.5*x*(j(n-1,x)-j(n+1,x) - j(n,x)/x);
}

dcomplex dhh(int n, dcomplex x){
  return h(n,x) + 0.5*x*(h(n-1,x)-h(n+1,x) - h(n,x)/x);
}


// *********** Polarizabilities *********************


dcomplex tE(int l, dfloat w, dfloat a){
dfloat x0 = w*a/Cspeed;
dcomplex xi = x0*sqrt(eps_dfloat(w));
return -imag_unit*( eps_dfloat(w)*j(l, xi)*djj(l,x0) - j(l,x0)*djj(l,xi) )/( h(l,x0)*djj(l,xi) - eps_dfloat(w)*j(l,xi)*dhh(l,x0) ); 
}

dcomplex tM(int l, dfloat w, dfloat a){
dfloat x0 = w*a/Cspeed;
dcomplex xi = x0*sqrt(eps_dfloat(w));
return -imag_unit*( j(l, xi)*djj(l,x0) - j(l,x0)*djj(l,xi) )/( h(l,x0)*djj(l,xi) - j(l,xi)*dhh(l,x0) );
}


dcomplex fs(int l, dcomplex z){
dfloat dl = dfloat(1)*l;
  return (dl + dfloat(1))*imag_unit*h(l,z)/z  - imag_unit*h(l+1,z);
}

dcomplex fe(int l, dcomplex z){
dfloat dl = dfloat(1)*l;
  return (dl + dfloat(1))*j(l,z)/z  - j(l+1,z);
}


//********************************************************************************
// ****************** Gauss Konrod Quadrature weights and abscisas ********************
//********************************************************************************

//***********************************************************************
void gauss_konrod_w(double GKQt[2*nw1 + 1][3]){

  int const N1 = nw1;

   auto XGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::abscissa();    
   auto WGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::weights();      
   auto WGt =  boost::math::quadrature::gauss<double, N1>::weights();

   double WGLt[2*N1+1];

// Changes the n order array of Gaussian quadrature to a 2n+1 array
for(int i = 0; i < XGKt.size(); ++i){
      if (i % 2 == 0){
           WGLt[i] = WGt[i/2];}
      else WGLt[i] = 0.;
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
}



void Omegas(dfloat xi[2*nw1 + 1], dfloat xk[2*nw1 + 1], dfloat xg[2*nw1 + 1]){

double GKQw1[2*nw1 + 1][3];
gauss_konrod_w(GKQw1);

double w;
int NN1 = 2*nw1 + 1;

for(int i=0; i < NN1;  ++i){
xi[i] = ((w2-w1)/2.)*GKQw1[i][0] + ((w2+w1)/2.);;
xk[i] = ((w2-w1)/2.)*GKQw1[i][1];
xg[i] = ((w2-w1)/2.)*GKQw1[i][2];
 }
}








void DL(int Lmax, dfloat r, dfloat vv, dfloat b, dfloat a, dcomplex DLy[6]){

dfloat dm;
dfloat dl;

dfloat dldwy[6];
double ddldwy[6];
double dvv = static_cast<double>(vv);
double dw =0.0;

dcomplex CM[Lmax][2*Lmax +1], tMl[Lmax];
dcomplex DE[Lmax][2*Lmax +1], tEl[Lmax];
dcomplex dz[2][Lmax], df[2][Lmax];

dcomplex IErty[4], IHrty[4];
dcomplex IErfy[4], IHrfy[4];

dfloat dl1, dl2, dm1, dm2;
double normfactor; // A normalization factor due to the normalization of integrals

double IN1, IV1, IW1, IW2, IU1, IU2;
FUN(IM,IN,IU,IV,IW,IX,IY,IZ,ID,III);


dcomplex AA[Lmax][2*Lmax+1], BB[Lmax][2*Lmax+1];
dcomplex aa, bb;


// Calls the Gauss - Konrod function
dfloat xi[2*nw1 + 1], xk[2*nw1 + 1], xg[2*nw1 + 1];
Omegas(xi, xk, xg);

dfloat w, k0;

char filenamex[sizeof "Results/50/quad_dldw_a1nm_v0.99c_b1.5nm_extF2.dat"];
sprintf(filenamex, "Results/%s/%d/quad_dldw_a%.2gnm_v%.2g_b%.2gnm.dat", directoryLv, Lmax, static_cast<double>(a/(1.*nm)), static_cast<double>(vv), static_cast<double>(b/(1.*nm)));
FILE *fpx = fopen(filenamex,"w+");
fprintf(fpx,"Angular Momentum Spectrum, a: %.2gnm    v: %.2gc   b: %.2gnm    Lmax: %d \n", static_cast<double>(a/(1.*nm)), static_cast<double>(vv), static_cast<double>(b/(1.*nm)));
fprintf(fpx,"\n");
fprintf(fpx,"         betav         w(au)                   dlEdwy                   dlHdwy                 dlEsdwy                   dlHsdwy                   dlEsdwyExt                   dlHsdwyExt\n");

for (int l = 1; l <= Lmax; ++l){
    for (int m = -l; m <= l; ++m){ 
    AA[l-1][l+m] = A(l,m,vv);
    BB[l-1][l+m] = B(l,m,vv);
    }
}


for (int i = 0; i < 6; ++i){
DLy[i] = dfloat(0);
}




for (int i = 0; i < iw1; ++i){     // Here we can select the Omega interval to integrate "2*NN + 4" for total.
    w = dfloat(1)*xi[i];
    dw= static_cast<double>(xi[i]);
    k0 = w/Cspeed;


// Calculates  common parameters of the scatterred fields
for (int l = 1; l <= Lmax; l++){
    tEl[l-1] = tE(l,w,a);
    tMl[l-1] = tM(l,w,a);
    
    dz[0][l-1] = imag_unit*h(l,k0*r);
    df[0][l-1] = fs(l,k0*r);
    dz[1][l-1] = j(l,k0*r);
    df[1][l-1] = fe(l,k0*r);

    for (int m = -l; m <= l; ++m){

        aa = AA[l-1][m+l];
        bb = BB[l-1][l+m];
    
        DE[l-1][m + l] = pow(imag_unit,l)*alm(l,m)*PsiE(l,m,w,b,vv, bb);
        CM[l-1][m + l] = pow(imag_unit,l)*alm(l,m)*PsiM(l,m,w,b,vv, aa);
    
    }
} 


// Calculates  common parameters of the scatterred fields
for (int l = 1; l <= Lmax; l++){
    tEl[l-1] = tE(l,w,a);
    tMl[l-1] = tM(l,w,a);
    
    dz[0][l-1] = imag_unit*h(l,k0*r);
    df[0][l-1] = fs(l,k0*r);
    dz[1][l-1] = j(l,k0*r);
    df[1][l-1] = fe(l,k0*r);

    for (int m = -l; m <= l; ++m){

        aa = AA[l-1][m+l];
        bb = BB[l-1][l+m];
    
        DE[l-1][m + l] = pow(imag_unit,l)*PsiE(l,m,w,b,vv, bb);
        CM[l-1][m + l] = pow(imag_unit,l)*PsiM(l,m,w,b,vv, aa);

        //DE[l-1][m + l] = pow(imag_unit,l)*alm(l,m)*PsiE(l,m,w,b,vv, bb);
        //CM[l-1][m + l] = pow(imag_unit,l)*alm(l,m)*PsiM(l,m,w,b,vv, aa);
    
    }
} 

for (int rr = 0; rr < 6; ++rr){
dldwy[rr] = dfloat(0);
ddldwy[rr] = 0.0;
}


for (int l1 = 1; l1 <= Lmax; ++l1){
     dl1 = dfloat(1)*l1;
     // Calculate progress percentage
                float progress = 100.0f * (l1 / float(Lmax));
                // Print dynamic line
                std::cout << "\r"
                          << "Progress for i=" << i <<" out of "<< iw1 <<" : " 
                          << std::fixed << std::setprecision(1) << progress << "% "
                          << std::flush;

for (int l2 = 1; l2 <= Lmax; ++l2){
         dl2 = dfloat(1)*l2; 
         for (int m1 = -l1; m1 <= l1; ++m1){
              dm1 = dfloat(1)*m1;
              for (int m2 = -l2; m2 <= l2; ++m2){
                   dm2 = dfloat(1)*m2;

                    if(m2 == m1+1 || m2 == m1-1){

                                    // Radial - Radial 
                        normfactor = pow((2.*static_cast<double>(l1)+1.)*(static_cast<double>(l1+m1)
                            +1.)/((2.*static_cast<double>(dl1)+3.)*(static_cast<double>(l1-m1)+1.)),0.5);

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

                     IErty[rr] = dfloat(0);
                     IHrty[rr] = dfloat(0);

                     IErfy[rr] = dfloat(0);
                     IHrfy[rr] = dfloat(0);
                    }
                   }

                    dldwy[0] +=  (dfloat(1)/(4*Pi))*pow(r,3)*(( IErty[2] - IErfy[2] ).real()
                                                    +( IErty[3] - IErfy[3] ).real());
                    dldwy[1] +=  (dfloat(1)/(4*Pi))*pow(r,3)*((IHrty[2] -IHrfy[2]).real()
                                                    +(IHrty[3] -IHrfy[3]).real());
                    dldwy[2] +=  (dfloat(1)/(4*Pi))*pow(r,3)*( IErty[0] - IErfy[0] ).real();
                    dldwy[3] +=  (dfloat(1)/(4*Pi))*pow(r,3)*(IHrty[0] -IHrfy[0]).real();
                    dldwy[4] +=  (dfloat(1)/(4*Pi))*pow(r,3)*( IErty[1] - IErfy[1] ).real();
                    dldwy[5] +=  (dfloat(1)/(4*Pi))*pow(r,3)*(IHrty[1] -IHrfy[1]).real();
            } // for m2
        }  // for m1

} // for l2
} //for l1

// Clear the line after completion
std::cout << "\r" << std::string(50, ' ') << "\r" << std::flush;


for (int rr = 0; rr < 6; ++rr){ 
  ddldwy[rr] = static_cast<double>(dldwy[rr]);
  DLy[rr] += (xk[i]*(dfloat(1) + imag_unit) - imag_unit*xg[i])*dldwy[rr];
}

// Here print the dldw's, for each, w, l
fprintf(fpx,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",static_cast<double>(vv),static_cast<double>(w),ddldwy[0],ddldwy[1],ddldwy[2], ddldwy[3], ddldwy[4],ddldwy[5]);

}// for w

std::cout << std::setprecision(std::numeric_limits<typename dcomplex::value_type>::digits10);

cout << "v : " << vv << endl;

//cout << "DLEy : " << DLy[0] << endl;
//cout << "DLHy : " << DLy[1] << endl;                // prints result
//cout << endl;

//cout << "DLEsy : " << DLy[2] << endl;
//cout << "DLHsy : " << DLy[3] << endl;                // prints result
//cout << endl;

cout << "DL : " << DLy[0] + DLy[1] + DLy[2] + DLy[3] << endl;
//cout << endl;

// Close the files
fclose(fpx);

} //end void





int main(void){


double nm1 = 100./5.2918;

double v;                    // Defines parameter electron speed vv, frequency a.u,
double b= 51;  
double a = 50;                        // impact parameter b, NP ratio, Integration ratio r.
double r = 50.05;

dcomplex DLy[6];
dfloat Ly;

int Lmax=51;

char filename[sizeof "Results/50/50/quadHead_DL_a1nm_v0.5c_vs_v.dat"];
sprintf(filename, "Results/%s/%d/quadHead_DL_a%.2gnm_b%.2gnm_vs_v.dat", directoryLv, Lmax, static_cast<double>(a), static_cast<double>(b));
FILE *fpp = fopen(filename,"a+");
fprintf(fpp,"Total angular momentum transfered, a: %.2gnm     v: %.2gc    Lmax: %d  \n", a/nm, v, Lmax);
fprintf(fpp,"\n");
fprintf(fpp,"     b        DPx     DPEesX     DPHesX     DPEssX     DPHssX     DPEeeX     DPHeeX\n");

char filenamee[sizeof "Results/50/50/quadHead_DL_a1nm_v0.5c_vs_v_error.dat"];
sprintf(filenamee, "Results/%s/%d/quadHead_DL_a%.2gnm_b%.2gnm_vs_v_error.dat", directoryLv, Lmax, static_cast<double>(a), static_cast<double>(b));
FILE *fppe = fopen(filenamee,"a+");
fprintf(fppe,"Total momentum transfered err X, a: %.2gnm     v: %.2gc    Lmax: %d  \n", a, b, Lmax);
fprintf(fppe,"\n");
fprintf(fppe,"     b        DPx     DPEesX     DPHesX     DPEssX     DPHssX     DPEeeX     DPHeeX\n");


int k, rr;
double ddLy, dly[4],idly[4];



double vv[] = {0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95} ;


cout <<  "Starting ... "<< endl << endl;

//#pragma omp parallel for shared(r, b, a,vv) private(b,k,rr,Px,Pz,DPx,DPz) num_threads(8)
for( k=0; k < 10; ++k){

v = vv[k];

DL(Lmax, r*nm, dfloat(1)*v, b*nm, a*nm, DLy); 

Ly=0.;

for (rr = 0; rr < 6; ++rr){ 
  Ly += DLy[rr].real();

  dly[rr] = static_cast<double>(DLy[rr].real());

  idly[rr] = static_cast<double>(DLy[rr].imag());

}


ddLy = static_cast<double>(Ly);

fprintf(fpp ,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",v,ddLy,    dly[0], dly[1], dly[2], dly[3], dly[4], dly[5]);
fprintf(fppe ,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",v,     idly[0],idly[1],idly[2],idly[3],idly[4],idly[5]);

}
            



return(0);
}
//**************************************************** END ****************************************************