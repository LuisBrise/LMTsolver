//*****************************************************************************************************************
//
//   DLGK_CALIBRATE.cpp: NOT TESTED
//
//    This program calculates the spectral angular momentum dl/dw transferred to a 
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
//#include "AMT_CALIBRATE_cut.h"

// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|

const int Lmax = 13;

//************************************************

double Cspeed = 137.035999139;
double Pi     = boost::math::constants::pi<double>();           // Parameter for Drude model of the NP and atomic units
double nm     = 100./5.2918;

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
	//res += pow(betav,-(l+1))*pow(1i,l-j)*factorial2(2*l+1)*alm(l,m)*II(l,m,j,l-j)
	///(pow(2.*gam,j)*factorial(l-j)*factorial((j-m)/2)*factorial((j+m)/2));
	res += pow(betav,-(l+1))*pow(1i,l-j)*factorial2(2*l+1)*II(l,m,j,l-j)
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


//*********************************************************************************************************
// Angular momentum spectrum:
//
// [0] Electric contribution. [1] Magnetic contribution.  Cross
// [2] Electric contribution. [3] Magnetic contribution.  Scat
// [4] Electric contribution. [5] Magnetic contribution.  Ext
//*********************************************************************************************************
//*********************************************************************************************************
void dldw(double r, double w, double vv, double b, double a, double dldwy[6]){

dcomplex CM[Lmax][2*Lmax +1], tMl[Lmax];
dcomplex DE[Lmax][2*Lmax +1], tEl[Lmax];
dcomplex dz[2][Lmax], df[2][Lmax];

dcomplex IErty[4], IHrty[4];
dcomplex IErfy[4], IHrfy[4];

double dl1, dl2, dm1, dm2;

double IN1, IV1, IW1, IW2, IU1, IU2;
FUN(IM,IN,IU,IV,IW,IX,IY,IZ,ID,III);


dcomplex AA[Lmax][2*Lmax+1], BB[Lmax][2*Lmax+1];
dcomplex aa, bb;

for (int l = 1; l <= Lmax; ++l){
    for (int m = -l; m <= l; ++m){ 
    AA[l-1][l+m] = A(l,m,vv);
    BB[l-1][l+m] = B(l,m,vv);
    }
}

double k0 = w/Cspeed;

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

                        IV1 = IV[l1-1][l2-1][m1+l1][m2+l2];

                        IW1 = IW[l1-1][l2-1][m1+l1][m2+l2];
                        IW2 = IW[l1][l2-1][m1+l1+1][m2+l2];   

                        IU1 = IU[l1-1][l2-1][m1+l1][m2+l2];
                        IU2 = IU[l1][l2-1][m1+l1+1][m2+l2];

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

}

double dl(double a, double r, double b, double vv, double x){

double dldwy[6];

if( BesselK(0,x) == 0. ){
  return 0.;}
else{
dldw(r, x, vv, b, a, dldwy);      // Calls function momentum dpdw
return (dldwy[0] + dldwy[1] + dldwy[2] + dldwy[3]);}
}

//*********************************************************************************************************
//********************************************** MAIN PROGRAM *********************************************
int main(void){

double vv = 0.5;
double a = 5.*nm;                        // impact parameter b, NP ratio, Integration ratio r.
double r = a + 0.05*nm;
double b = a + 0.5*nm;

double lf = 260.;  // last freq ---> 40, v0.3 b1.5;  65, v0.5 b1.5;  100, v0.3 b1.5;  160, v0.9 b1.5; 

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