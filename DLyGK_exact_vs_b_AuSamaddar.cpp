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
//      g++ -o DLyGK_exact_vs_b.out DLyGK_exact_vs_b.cpp -lcomplex_bessel -w 
//      g++ -o DLyGK_exact_vs_b_Au.out DLyGK_exact_vs_b_Au.cpp -lcomplex_bessel -w && ./DLyGK_exact_vs_b_Au.out 
//
//          By Jesús Castrejon, jcastrejon@ciencias.unam.mx (25/02/2019)
// Modified by Jorge Luis Briseño, jorgeluisbrisenio@ciencias.unam.mx (21/02/2023)
//
//*****************************************************************************************************************
//*****************************************************************************************************************


#include "IN31.h"

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
double wp     = 0.48297106037599996;
double Gamma  = 0.007244389509;

//************************************************
// Maximum number of multipole to take into account

int Lmax = 3;

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
                   // In Au, b=1.5nm a=1nm, for v=0.5 appx 35,  v=0.8 appx    ,v=.99 appx 180..

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

double IN1, IV1, IW1, IW2, IU1, IU2;
double IK1, IK2, IK3, IK4, nnpp, nnp, nnm, nnmm, dlnpp, dlnp, dlnm, dll, dllp;
FUN(IM,IN,IU,IV,IW,IX,IY,IZ,ID,III);


// Calls the Gauss - Konrod function
double xi[2*NN + 4], xk[2*NN + 4], xg[2*NN + 4];
Omegas(xi, xk, xg);

double w, k0;

char filenamex[sizeof "Results/dldw_a1nm_v0.99c_b1.5nm_extF2.dat"];
sprintf(filenamex, "Results/LybAu/%d/dldw_a%.2gnm_v%.2g_b%.2gnm.dat", Lmax, a/(1.*nm), vv, b/(1.*nm));
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
    
        DE[l-1][m + l] = pow(1i,l)*alm(l,m)*PsiE(l,m,w,b,vv, bb);
        CM[l-1][m + l] = pow(1i,l)*alm(l,m)*PsiM(l,m,w,b,vv, aa);
    
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

                        if(m2 == m1 +1){
                           if(m1 < 0){
                              modf( (dl1-dl2+1.0)/2 , &nnp );
                              modf( (dl1-dl2-2.0)/2 , &nnm );
                              modf( (dl1-dl2-1.0)/2 , &nnpp );
                              if(nnp<0){nnp=0.0;}
                              if(nnm<0){nnm=0.0;}
                              if(nnpp<1){nnpp=-1.0;}
                              dlnp=0.0; // delta_{l2, l1-2*nnp + 1}
                              dlnm=0.0; //delta_{l2, l1-2*nnm-2}
                              dll=0.0; // delta_{l1,l2}
                              dllp=0.0; // delta_{l1,l2+1}
                              dlnpp=0.0; // delta_{l1,l2+2nnpp+1}
                              if(dl2==dl1-2.*nnp+1.){dlnp=1.0;}
                              if(dl2==dl1-2.*nnm-2.){dlnm=1.0;}
                              if(dl2==dl1){dll=1.0;}
                              if(dl2==dl1-1){dllp=1.0;}
                              if(dl2==dl1-2*nnpp-1){dlnpp=1.0;}
                              IK1 = 2.0*factorial(l1+m1)*dlnp/factorial(l1-m1);
                              IK2 = 2.*m1*factorial(l1+m1)*dlnm/factorial(l1-m1)
                                    -2.*l1*factorial(l2+m2)*dll/( (2.*l1+1)*factorial(l1-m1) );
                              IK3 = 2.*factorial(l1+m1)*dlnm/factorial(l1-m1)
                                    +2.*factorial(l2+m2)*dll/( (2.*l1+1)*factorial(l1-m1) );
                              IK4 = -factorial(l1+m1)*( ( (l2+m2-1)(l2+m2)-(l2-m2)(l2-m2+1) )*dlnpp
                                    + ( (l1+m1)(l1+m1-1) + (l1-m1)(l1-m1-1)/(2.*l1+1) )*dllp )
                                    /( (2.*l2+1)*factorial(l1-m1) );
                           }
                           else{
                              dll = 0.0;
                              dlnp = 0.0;
                              if(l1 == l2){dll=1.0;}
                              if(l1 == l2 + 1){dlnp = 1.0;}
                              // (2.*l1+1.)*factorial(l1-m1)/(4.*Pi*factorial(l1+m1))
                              IK1 = 0;
                              IK2 = -2.0*(l1+1.)*factorial(l1+m1)*dll/( (2.*l1+1.)*factorial(l1-m1-1) );
                              IK3 = -2.0*factorial(l1+m1)*dll/( (2.*l1+1.)*factorial(l1-m1-1) );
                              //pow((2.*l1+1.)*factorial(l2-m2)/(4.*Pi*factorial(l2+m2+2)),0.5)*pow((2.*l2+1.)*factorial(l2-m2)/(4.*Pi*factorial(l2+m2)),0.5)
                              IK4 = -2.0*(l1+1.)*factorial(l2+m2)*dlnp/( (2.*l1+1)*(2.*l2+1)*factorial(l2-m2) );
                           }
                        }
                        else{
                           if(m1<0){
                              dll=0.0;
                              dlnp=0.0;
                              if(l1 == l2){dll=1.0;}
                              if(l1 == l2 + 1){dlnp = 1.0;}
                              IK1=0;
                              IK2=2.0*factorial(l1+m1)*(l1+1)*dll/ ( factorial(l1-m1+1.)*(2.*l1+1.) );
                              IK3=2.0*factorial(l1+m1)*dll/( factorial(l1-m1+1)*(2.*l1+1.) );
                           }
                           else{
                              modf( (dl1-dl2-1.0)/2, &nnm )
                              if(nnm < 0){
                                 nnm=0; // only positive integers are allowed in Samaddar solution
                              }
                              dlnm=0.0;
                              if(dl2 == dl1 - 2.0*nnm - 1.0){
                                 dlnm=1.0;
                              }
                              IK1 = -2.0*factorial(l2+m1)*dlnm/factorial(l2-m2)
                           }
                        } 

                        // IV kernel x^2/sqrt(1-x^2)
                        IV1 = IV[l1-1][l2-1][m1+l1][m2+l2];

                        // IW kernel x/sqrt(1-x^2)
                        IW1 = IW[l1-1][l2-1][m1+l1][m2+l2];
                        IW2 = IW[l1][l2-1][m1+l1+1][m2+l2];   

                        // IU kernel 1/sqrt(1-x^2)
                        IU1 = IU[l1-1][l2-1][m1+l1][m2+l2];
                        IU2 = IU[l1][l2-1][m1+l1+1][m2+l2];

                        /* Integrants, 0 -> Scat
                                       1 -> Ext
                                       2 -> Ext-Scat 
                                       3 -> Scat-Ext (second is conjugated)*/



                    // Zenital - Radial

                    // y component

                     // pow((2.*l1+1.)*factorial(l1-m1)/(4.*Pi*factorial(l1+m1)),0.5)*pow((2.*l2+1.)*factorial(l2-m2)/(4.*Pi*factorial(l2+m2)),0.5)
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

cout << "b : " << b/nm << endl;

cout << "DLEy : " << DLy[0] << endl;
cout << "DLHy : " << DLy[1] << endl;                // prints result
cout << endl;

cout << "DLEsy : " << DLy[2] << endl;
cout << "DLHsy : " << DLy[3] << endl;                // prints result
cout << endl;

cout << "DL : " << DLy[0] + DLy[1] + DLy[2] + DLy[3] << endl;
cout << endl; 

} //end void



//*********************************************************************************************************
//********************************************** MAIN PROGRAM *********************************************
int main(void){


cout.precision(17);

cout << endl;
cout << "Momentum Gauss - Kronrod :" << endl;
cout << endl;
//cout << "Lmax = " << Lmax << endl;
//cout << endl;

/*double b = 21.5*nm;                        
double a = 20.0*nm;                        
double r = 20.05*nm;
double vv = 0.5;*/

/*double b = 51.5*nm;                        
double a = 50.0*nm;                        
double r = 50.05*nm;
double vv = 0.5;*/

double b = 5.5*nm;                        
double a = 5.0*nm;                        
double r = 5.05*nm;
double vv = 0.5;


/*double b = 11.5*nm;                        
double a = 10.0*nm;                        
double r = 10.05*nm;
double vv = 0.5;*/

dcomplex DLy[6];
double Ly;

// Start the timer
auto start = std::chrono::high_resolution_clock::now();
for (Lmax = 1; Lmax < 30; Lmax++){
cout << "Lmax = " << Lmax << endl;
cout << endl;
char filename[sizeof "Results/DL_a1nm_v0.99c_b1.5nm_extF2.dat"];
sprintf(filename, "Results/LybAu/%d/DL_a%.2gnm_v%.2g.dat", Lmax,a/(1.*nm), vv);
FILE *fpp = fopen(filename,"w+");
fprintf(fpp,"Total Angular momentum transfered, a: %.2gnm    v: %.2gc   Lmax: %d  \n", a/(1.*nm), vv, b/(1.*nm), Lmax);
fprintf(fpp,"\n");
fprintf(fpp,"         b         DLEy                  DLHy                  DLEsy                  DLHsy                  DLEsyExt                  DLHsyExt                  DLy\n");


char filenamer[sizeof "Results/DL_a1nm_v0.99c_b1.5nm_error_extF2.dat"];
sprintf(filenamer, "Results/LybAu/%d/DL_a%.2gnm_v%.2g_error.dat", Lmax, a/(1.*nm), vv);
FILE *fppe = fopen(filenamer,"w+");
fprintf(fppe,"Total Angular momentum transfered (error), a: %.2gnm    v: %.2gc   Lmax: %d \n", a/(1.*nm), vv, b/(1.*nm), Lmax);
fprintf(fppe,"\n");
fprintf(fppe,"         b       errDLEy                errDLHy                 errDLEsy                errDLHsy                 errDLEsyExt                errDLHsyExt\n");

for (int i = 0; i <= 9; ++i)
{
   b = (5.5 + (10.0-5.5)*i/9)*nm;
   //b = (11.5 + (20.0-11.5)*i/10)*nm;
   //b = (21.5 + (30.0-21.5)*i/10)*nm;
   //b = (51.5 + (60.0-51.5)*i/10)*nm;
   DL(r, vv, b, a, DLy);
   Ly=0.0;

   for (int rr = 0; rr < 6; ++rr){ 
   Ly += DLy[rr].real();}

   // Here print the total momentum
   fprintf(fpp,"%.17g %.17g %.17g %.17g %.17g  %.17g %.17g %.17g \n",b/(1.*nm),DLy[0].real(),DLy[1].real(),DLy[2].real(),DLy[3].real(),DLy[4].real(),DLy[5].real(), Ly);
   fprintf(fppe,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",b/(1.*nm),DLy[0].imag(),DLy[1].imag(),DLy[2].imag(),DLy[3].imag(),DLy[4].imag(),DLy[5].imag());
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
