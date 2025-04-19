#include "CoreTypes.h"
#include "CoreIncludes.h"
#include "BoostIncludes.h"
#include "IN51.h"                

//**********************************************************************
// Dielectric Function
//********************************************************************** 
#include "Epsilon.h"

// IM, IN, IU, IV, IW, IX, IY, IZ, ID, II
static Irreducible_Integrals integrals;

//This is where we are gonna input Gauss abscissas and weights
double xi[2*NN + 4], xk[2*NN + 4], xg[2*NN + 4];

#include "CommonFunctions.h"

#include "Quadrature.h"

#include "GeneratePaths.h"

void print_initial_message(SimulationConfig& config){
    cout << 
R"(LMTsolver

Description:
This program calculates the spectral linear momentum (dP/dω) transferred to a nanoparticle (NP) due to the interaction with a fast-moving electron. The calculation is performed by analytically solving the closed surface integral, resulting in a double multipolar sum expression. The program focuses exclusively on the external field contribution.This program calculates the total angular momentum transfer (DeltaP) for different values of the impact parameter b.

Key Features:
- The frequency integral is computed using the Gauss-Kronrod quadrature method, with an adaptive partition for rapid convergence.
- The spectral results are saved to files named "dpdw_*.dat".
- The integrated angular momentum results are saved to files named "DP_*.dat".

Usage:
- The user can input parameters through the terminal.
- The program supports multiple parameter sets and stops when the solution converges to a specified relative error threshold or reaches the maximum multipolar order (Lmax = 50).

Parameters:)" << endl;
cout << "Material : " << config.material << endl
     << "Lmax : "   << config.Lmax << endl             
     << "NP radius (nm) : " << config.a << endl      
     << "Impact parameter (nm) : " << config.b << endl
     << "Speed (c) : " << config.velocity << "c" << endl
     << "Is it a V Scan? - " << config.isVScan << endl
     << "Is it a B Scan? - " << config.isBScan << endl
     << "Is it a V vs B contour? - " << config.isBvsVContour << endl
     << "Multipolar threshold : " << config.errorThreshold << endl
     << "Number of points in the study : " << config.numOfPoints << endl
     << "Minimum value of Lmax : " << config.minLmax << endl
     << "______________________________________________________" << endl;
}

void DP(SimulationConfig& config, dcomplex DPx[6], dcomplex DPz[6]){

double r = config.r*nm;
double vv = config.velocity;
double b = config.b*nm;
double a = config.a*nm;
int Lmax = config.Lmax;

double dm;
double dl;

double dpdwx[6], dpdwz[6];

dcomplex CM[Lmax][2*Lmax +1], tMl[Lmax];
dcomplex DE[Lmax][2*Lmax +1], tEl[Lmax];
dcomplex dz[2][Lmax], df[2][Lmax];

dcomplex IErrx[4], IErrz[4], IHrrx[4], IHrrz[4];
dcomplex IErtx[4], IErtz[4], IHrtx[4], IHrtz[4];
dcomplex IErf[4], IHrf[4];
dcomplex IEttx[4], IEttz[4], IHttx[4], IHttz[4];
dcomplex IEffx[4], IEffz[4], IHffx[4], IHffz[4];

double dl1, dl2, dm1, dm2;
double normfactor2,normfactor3,normfactor4; // A normalization factor due to the normalization of integrals

double IN1, IV1, IW1, IW2, IW3, IU1, IU2, IU3, IU4;
double IM1, IZ1, IX1, IX2, IX3, IY1, IY2, IY3, IY4, ID1, ID2;

double w, k0;

auto filenamex = generate_dpdw_path(config, "x");
std::ofstream outx(filenamex);
outx.precision(17);
outx << "vv\t\t\tb\t\t\tw(au)\t\t\tdpEdwxInt\t\t\tdpHdwxInt\t\t\tdpEdwxScat\t\t\tdpHdwxScat\t\t\tdpEdwxExt\t\t\tdpHdwxExt\n";

auto filenamez = generate_dpdw_path(config, "z");
std::ofstream outz(filenamez);
outz.precision(17);
outz << "vv\t\t\tb\t\t\tw(au)\t\t\tdpEdwzInt\t\t\tdpHdwzInt\t\t\tdpEdwzScat\t\t\tdpHdwzScat\t\t\tdpEdwzExt\t\t\tdpHdwzExt\n";

dcomplex aa, bb;

for (int i = 0; i < 6; ++i){
DPx[i] = 0.;
DPz[i] = 0.;
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
    
    }
} 



for (int rr = 0; rr < 6; ++rr){
dpdwx[rr] = 0.;
dpdwz[rr] = 0.;
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
                        normfactor2 = pow((2.*dl1+1.)*(l1+m1+1.)/((2.*dl1+3.)*(l1-m1+1.)),0.5);
                        normfactor3 = pow((2.*dl2+1.)*(l2+m2+1.)/((2.*dl2+3.)*(l2-m2+1.)),0.5);
                        normfactor4 = normfactor2*normfactor3;

                        IN1 = integrals.IN[l1-1][l2-1][m1+l1][m2+l2];
                        IV1 = integrals.IV[l1-1][l2-1][m1+l1][m2+l2];

                        IW1 = integrals.IW[l1-1][l2-1][m1+l1][m2+l2];
                        IW2 = normfactor2*integrals.IW[l1][l2-1][m1+l1+1][m2+l2];   
                        IW3 = normfactor3*integrals.IW[l1-1][l2][m1+l1][m2+l2+1];

                        IU1 = integrals.IU[l1-1][l2-1][m1+l1][m2+l2];
                        IU2 = normfactor2*integrals.IU[l1][l2-1][m1+l1+1][m2+l2];
                        IU3 = normfactor3*integrals.IU[l1-1][l2][m1+l1][m2+l2+1];
                        IU4 = normfactor4*integrals.IU[l1][l2][m1+l1+1][m2+l2+1];

                        /* Integrants, 0 -> Scat
                                       1 -> Ext
                                       2 -> Ext-Scat 
                                       3 -> Scat-Ext (second is conjugated)*/


                    IErrx[0] = tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; 

                    IHrrx[0] = tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; 


                    IErrx[1] = DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; 

                    IHrrx[1] = CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; 



                    IErrx[2] = DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; //G

                    IHrrx[2] = CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; //G



                    IErrx[3] = tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; //G

                    IHrrx[3] = tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; //G


                                   // Zenital - Radial

                    IErtx[0] = -(tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 - (dl1-dm1+1.)*IW2)
                               + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);
     
                    IHrtx[0] =  (- tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1- (dl1-dm1+1.)*IW2)
                              +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);



                    IErtx[1] = -(DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 - (dl1-dm1+1.)*IW2)
                               + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);
     
                    IHrtx[1] =  (- CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1- (dl1-dm1+1.)*IW2)
                              +    DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);



                    IErtx[2] = -(DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 - (dl1-dm1+1.)*IW2)
                               + CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);
     
                    IHrtx[2] =  (- CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1- (dl1-dm1+1.)*IW2)
                              +    DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);




                    IErtx[3] = -(tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 - (dl1-dm1+1.)*IW2)
                               + tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);
     
                    IHrtx[3] =  (- tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1- (dl1-dm1+1.)*IW2)
                              +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);




                                    // Azimutal - Radial


                    IErf[0] =  (dm2-dm1)*( tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                      +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1); 

                    IHrf[0] =  (dm2-dm1)*( -tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                          + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);  

                    IErf[1] =  (dm2-dm1)*( CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                      +    DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1); 

                    IHrf[1] =  (dm2-dm1)*( -DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                          + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);  

                    IErf[2] =  (dm2-dm1)*( CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                      +    DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1); 

                    IHrf[2] =  (dm2-dm1)*( -DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                          + CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);

                    IErf[3] =  (dm2-dm1)*( tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                      +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1); 

                    IHrf[3] =  (dm2-dm1)*( -tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                          + tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);


                                        // Azimutal - Azimutal

                    IEffx[0] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*dm1*dm2*IU1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); //G

                    IHffx[0] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*dm1*dm2*IU1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);   //   

                    IEffx[1] =  DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*dm1*dm2*IU1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); //G

                    IHffx[1] =  CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*dm1*dm2*IU1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);   //  



                    IEffx[2] =  DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*dm1*dm2*IU1 
                              + CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); //G

                    IHffx[2] =  CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*dm1*dm2*IU1 
                              - DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              - CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);   //  

                    IEffx[3] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*dm1*dm2*IU1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); //G

                    IHffx[3] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*dm1*dm2*IU1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);   //  




                                        //  Zenital - Zenital

                    IEttx[0] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*dm1*dm2*IU1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 

                    IHttx[0] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*dm1*dm2*IU1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 



                    IEttx[1] =  CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*dm1*dm2*IU1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 

                    IHttx[1] =  DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*dm1*dm2*IU1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);  


                    IEttx[2] =  CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*dm1*dm2*IU1 
                              + CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 

                    IHttx[2] =  DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*dm1*dm2*IU1 
                              - DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              - CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 


                    IEttx[3] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*dm1*dm2*IU1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 

                    IHttx[3] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*dm1*dm2*IU1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); }


                else{
                  for (int rr = 0; rr < 4; ++rr){
                     IErrx[rr] = 0.;
                     IHrrx[rr] = 0.;
                     IErtx[rr] = 0.;
                     IHrtx[rr] = 0.;
                     IEffx[rr] = 0.;
                     IHffx[rr] = 0.;
                     IErf[rr] = 0.;
                     IHrf[rr] = 0.;
                     IEttx[rr] = 0.;
                     IHttx[rr] = 0.;
                    }
                   }




                if(m2 == m1){
                                              // Radial - Radial 
                        normfactor2 = pow((2.*dl1+1.)*(l1+m1+1.)/((2.*dl1+3.)*(l1-m1+1.)),0.5);
                        normfactor3 = pow((2.*dl2+1.)*(l2+m2+1.)/((2.*dl2+3.)*(l2-m2+1.)),0.5);
                        normfactor4 = normfactor2*normfactor3;

                        IM1 = integrals.IM[l1-1][l2-1][m1+l1][m2+l2];
                        IZ1 = integrals.IZ[l1-1][l2-1][m1+l1][m2+l2];

                        ID1 = integrals.ID[l1-1][l2-1][m1+l1][m2+l2];
                        ID2 = normfactor2*integrals.ID[l1][l2-1][m1+l1+1][m2+l2];

                        IX1 = integrals.IX[l1-1][l2-1][m1+l1][m2+l2];
                        IX2 = normfactor2*integrals.IX[l1][l2-1][m1+l1+1][m2+l2];
                        IX3 = normfactor3*integrals.IX[l1-1][l2][m1+l1][m2+l2+1];

                        IY1 = integrals.IY[l1-1][l2-1][m1+l1][m2+l2];
                        IY2 = normfactor2*integrals.IY[l1][l2-1][m1+l1+1][m2+l2];
                        IY3 = normfactor3*integrals.IY[l1-1][l2][m1+l1][m2+l2+1];
                        IY4 = normfactor4*integrals.IY[l1][l2][m1+l1+1][m2+l2+1];



                    IErrz[0] = tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; // G of "Good" cheked with TestGKL.cpp

                    IHrrz[0] = tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; //G


                    IErrz[1] = DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; // G of "Good" cheked with TestGKL.cpp

                    IHrrz[1] = CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; //G


                    IErrz[2] = DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; // G of "Good" cheked with TestGKL.cpp

                    IHrrz[2] = CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; //G


                    IErrz[3] = tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; // G of "Good" cheked with TestGKL.cpp

                    IHrrz[3] = tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; //G


                                              // Zenital - Radial

                    IErtz[0] = -( tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1);  //G
   
                    IHrtz[0] = (- tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1); //G

                    IErtz[1] = -( DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1);  //G
   
                    IHrtz[1] = (- CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1); //G

                    IErtz[2] = -( DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1);  //G
   
                    IHrtz[2] = (- CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1); //G

                    IErtz[3] = -( tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1);  //G
   
                    IHrtz[3] = (- tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1); //G

                                           // Azimutal - Azimutal

                    IEffz[0] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*dm1*dm2*IY1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G

                    IHffz[0] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*dm1*dm2*IY1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G




                    IEffz[1] =  DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*dm1*dm2*IY1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G

                    IHffz[1] =  CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*dm1*dm2*IY1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G




                    IEffz[2] =  DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*dm1*dm2*IY1 
                              + CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G

                    IHffz[2] =  CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*dm1*dm2*IY1 
                              - DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              - CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G





                    IEffz[3] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*dm1*dm2*IY1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G

                    IHffz[3] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*dm1*dm2*IY1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G



                                                                      //  Zenital - Zenital

                    IEttz[0] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*dm1*dm2*IY1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); 

                    IHttz[0] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*dm1*dm2*IY1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); 




                    IEttz[1] =  CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*dm1*dm2*IY1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); 

                    IHttz[1] =  DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*dm1*dm2*IY1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  





                    IEttz[2] =  CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*dm1*dm2*IY1 
                              + CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); 

                    IHttz[2] =  DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*dm1*dm2*IY1 
                              - DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              - CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  


                    IEttz[3] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*dm1*dm2*IY1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); 

                    IHttz[3] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*dm1*dm2*IY1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); }  
                                              
                else{
                  for (int rr = 0; rr < 4; ++rr){
                     IErrz[rr] = 0.;
                     IHrrz[rr] = 0.;
                     IErtz[rr] = 0.;
                     IHrtz[rr] = 0.;
                     IEttz[rr] = 0.;
                     IHttz[rr] = 0.;
                     IEffz[rr] = 0.;
                     IHffz[rr] = 0.;
                    }
                   }




                     
                    dpdwx[0] += (1./(4.*Pi))*pow(r,2)*((0.5*(IErrx[2] - IEttx[2] - IEffx[2]) + IErtx[2] - IErf[2]).real()
                                            + (0.5*(IErrx[3] - IEttx[3] - IEffx[3]) + IErtx[3] - IErf[3])).real();
                    dpdwz[0] += (1./(2.*Pi))*pow(r,2)*((0.5*(IErrz[2] - IEttz[2] - IEffz[2]) - IErtz[2]).real()
                                            + (0.5*(IErrz[3] - IEttz[3] - IEffz[3]) - IErtz[3])).real();

                    dpdwx[1] +=  (1./(4.*Pi))*pow(r,2)*((0.5*(IHrrx[2] - IHttx[2] - IHffx[2]) + IHrtx[2] - IHrf[2]).real()
                                             + (0.5*(IHrrx[3] - IHttx[3] - IHffx[3]) + IHrtx[3] - IHrf[3])).real();
                    dpdwz[1] +=  (1./(2.*Pi))*pow(r,2)*((0.5*(IHrrz[2] - IHttz[2] - IHffz[2]) - IHrtz[2]).real()
                                             + (0.5*(IHrrz[3] - IHttz[3] - IHffz[3]) - IHrtz[3])).real();

                    dpdwx[2] +=  (1./(4.*Pi))*pow(r,2)*((0.5*(IErrx[0] - IEttx[0] - IEffx[0]) + IErtx[0] - IErf[0])).real();
                    dpdwz[2] +=  (1./(2.*Pi))*pow(r,2)*((0.5*(IErrz[0] - IEttz[0] - IEffz[0]) - IErtz[0])).real();

                    dpdwx[3] +=  (1./(4.*Pi))*pow(r,2)*((0.5*(IHrrx[0] - IHttx[0] - IHffx[0]) + IHrtx[0] - IHrf[0])).real();
                    dpdwz[3] +=  (1./(2.*Pi))*pow(r,2)*((0.5*(IHrrz[0] - IHttz[0] - IHffz[0]) - IHrtz[0])).real();

                    dpdwx[4] +=  (1./(4.*Pi))*pow(r,2)*((0.5*(IErrx[1] - IEttx[1] - IEffx[1]) + IErtx[1] - IErf[1])).real();
                    dpdwz[4] +=  (1./(2.*Pi))*pow(r,2)*((0.5*(IErrz[1] - IEttz[1] - IEffz[1]) - IErtz[1])).real();

                    dpdwx[5] +=  (1./(4.*Pi))*pow(r,2)*((0.5*(IHrrx[1] - IHttx[1] - IHffx[1]) + IHrtx[1] - IHrf[1])).real();
                    dpdwz[5] +=  (1./(2.*Pi))*pow(r,2)*((0.5*(IHrrz[1] - IHttz[1] - IHffz[1]) - IHrtz[1])).real();





            } // for m2
        }  // for m1

} // for l2
} //for l1

// Here print the dldw's, for each, w, l
outx << vv << '\t'
     << b << '\t'
     << dpdwx[0] << '\t'
     << dpdwx[1] << '\t'
     << dpdwx[2] << '\t'
     << dpdwx[3] << '\t'
     << dpdwx[4] << '\t'
     << dpdwx[5] << '\n';
outz << vv << '\t'
     << b << '\t'
     << dpdwz[0] << '\t'
     << dpdwz[1] << '\t'
     << dpdwz[2] << '\t'
     << dpdwz[3] << '\t'
     << dpdwz[4] << '\t'
     << dpdwz[5] << '\n';

for (int rr = 0; rr < 6; ++rr){ 
  DPx[rr] += (xk[i]*(1. + 1i) - 1i*xg[i])*dpdwx[rr];
  DPz[rr] += (xk[i]*(1. + 1i) - 1i*xg[i])*dpdwz[rr];
}



}// for w

// Close the files
outx.close();
outz.close();

//cout << "b : " << b/nm << " nm" <<endl;
cout << "DPx : " << (DPx[0] + DPx[1] + DPx[2] + DPx[3] + DPx[4] + DPx[5]).real() 
     << " +- " << (DPx[0] + DPx[1] + DPx[2] + DPx[3] + DPx[4] + DPx[5]).imag() << endl;
cout << "DPz : " << (DPz[0] + DPz[1] + DPz[2] + DPz[3] + DPz[4] + DPz[5]).real() 
     << " +- " << (DPz[0] + DPz[1] + DPz[2] + DPz[3] + DPz[4] + DPz[5]).imag() << endl; 

} //end void



void LMTsolver(SimulationConfig& config){

    print_initial_message(config);
    
    dcomplex DPx[6]; 
    dcomplex DPz[6];

    if (config.isVScan){

        // Vector to store Ly values for each i in the previous Lmax iteration
        vector<double> previousPx(config.numOfPoints + 1, 0.0); // Initialize to 0.0
        vector<double> previousPz(config.numOfPoints + 1, 0.0); // Initialize to 0.0

        auto multipolarfilename = generate_multconvan_path(config);
        std::ofstream mout(multipolarfilename);
        mout.precision(17);
        mout << "Lmax\t\t\tAverage Relative Error\n";

        // Start the timer
        auto start = std::chrono::high_resolution_clock::now();
        //
        // Start loop for Lmax
        //
        for (config.Lmax = config.minLmax; config.Lmax <= LSmax ; config.Lmax++){
            // Print dynamic line
            cout << "Progress: " 
                      << "Lmax : " << config.Lmax
                      << " | b : " << config.b << " nm"
                      << " | NP radius: " << config.a << " nm" << endl;
                      //
            auto filenamex = generate_LMT_path(config, "x");
            std::ofstream outx(filenamex);
            outx.precision(17);
            outx << "vv\tb\tDPExInt\t\t\tDPHxInt\t\t\tDPExScat\t\t\tDPHxScat\t\t\tDPExExt\t\t\tDPHxExt\t\t\tDPxTot\n";
            
            auto errorfilenamex = generate_LMT_path(config, "x", "error_");
            std::ofstream eoutx(errorfilenamex);
            eoutx.precision(17);
            eoutx << "vv\tb\terrDPExInt\t\t\terrDPHxInt\t\t\terrDPExScat\t\t\terrDPHxScat\t\t\terrDPExExt\t\t\terrDPHxExt\n";

            auto filenamez = generate_LMT_path(config, "z");
            std::ofstream outz(filenamez);
            outz.precision(17);
            outz << "vv\tb\tDPEzInt\t\t\tDPHzInt\t\t\tDPEzScat\t\t\tDPHzScat\t\t\tDPEzExt\t\t\tDPHzExt\t\t\tDPzTot\n";
            
            auto errorfilenamez = generate_LMT_path(config, "z", "error_");
            std::ofstream eoutz(errorfilenamez);
            eoutz.precision(17);
            eoutz << "vv\tb\terrDPEzInt\t\t\terrDPHzInt\t\t\terrDPEzScat\t\t\terrDPHzScat\t\t\terrDPEzExt\t\t\terrDPHzExt\n";

            double totalRelativeError = 0.0; // Accumulator for relative errors
            bool thresholdReached = false;

            for (int i = 0; i <= config.numOfPoints; ++i)
            {
                config.velocity = (config.vvInit 
                    + (config.vvFin-config.vvInit)*i/config.numOfPoints);
                //Initialize Scattering Functions AA and BB when vv once vv stops changing for LSmax
                Initialize_ScatteringFunctions(AA, BB, config.velocity);
                printf("v : %3.2f   -->   ", config.velocity);
                DP(config, DPx, DPz);
                
                double Px=0.0;
                double Pz=0.0;

                for (int rr = 0; rr < 6; ++rr){ 
                    Px += DPx[rr].real();
                    Pz += DPz[rr].real();}

                // Calculate the error
                double errorPx = fabs(Px - previousPx[i]);
                double errorPz = fabs(Pz - previousPz[i]);
                double relativeError = (0.5*errorPx / fabs(Px))+ (0.5*errorPz / fabs(Pz));

                // Accumulate the relative error
                totalRelativeError += relativeError;

               // Here print the total momentum
               // Here print the total momentum
                outx << config.velocity << '\t'
                     << config.b << '\t'
                     << DPx[0].real() << '\t'
                     << DPx[1].real() << '\t'
                     << DPx[2].real() << '\t'
                     << DPx[3].real() << '\t'
                     << DPx[4].real() << '\t'
                     << DPx[5].real() << '\t'
                     << Px << '\n';
                eoutx << config.velocity << '\t'
                      << config.b << '\t'
                      << DPx[0].imag() << '\t'
                      << DPx[1].imag() << '\t'
                      << DPx[2].imag() << '\t'
                      << DPx[3].imag() << '\t'
                      << DPx[4].imag() << '\t'
                      << DPx[5].imag() << '\n';
                outz << config.velocity << '\t'
                     << config.b << '\t'
                     << DPz[0].real() << '\t'
                     << DPz[1].real() << '\t'
                     << DPz[2].real() << '\t'
                     << DPz[3].real() << '\t'
                     << DPz[4].real() << '\t'
                     << DPz[5].real() << '\t'
                     << Pz << '\n';
                eoutz << config.velocity << '\t'
                      << config.b << '\t'
                      << DPz[0].imag() << '\t'
                      << DPz[1].imag() << '\t'
                      << DPz[2].imag() << '\t'
                      << DPz[3].imag() << '\t'
                      << DPz[4].imag() << '\t'
                      << DPz[5].imag() << '\n';
                // Update previousLy for the next iteration
                previousPx[i] = Px;
                previousPz[i] = Pz;
            }
            // Calculate the average relative error
            double averageRelativeError = totalRelativeError / (config.numOfPoints + 1);
            mout << config.Lmax << '\t' << averageRelativeError << '\n';
            cout << "Average Relative Error is: " << averageRelativeError << endl;
            cout << "------------------------------------------------------" << endl << endl;
            // Check if the average error is below the threshold
            if (averageRelativeError < config.errorThreshold) {
                thresholdReached = true;
                mout << "Convergence achieved at Lmax = " << config.Lmax
                         << " with average relative error = " << averageRelativeError << '\n';
                cout << "Convergence achieved at Lmax = " << config.Lmax 
                         << " with average relative error = " << averageRelativeError << endl;
                cout << "------------------------------------------------------" << endl << endl;
            }

            // Close the files
            outx.close();
            eoutx.close();
            outz.close();
            eoutz.close();

            // Check if the threshold was reached
            if (thresholdReached) {
                break; // Exit the outer loop
            }
                    // Clear the line after completion
            //std::cout << "\r" << std::string(50, ' ') << "\r" << std::flush;
        }
        // End the timer
        auto end = std::chrono::high_resolution_clock::now();
        // Calculate the duration
        std::chrono::duration<double> duration = end - start;
        // Print the duration in seconds
        cout << "Execution time: " << duration.count() << " seconds." << std::endl;
        mout << "Execution time: " << duration.count() << " seconds." << std::endl;
        mout.close();

    }

    else if (config.isBScan)
    {
        // Vector to store P values for each i in the previous Lmax iteration
        vector<double> previousPx(config.numOfPoints + 1, 0.0); // Initialize to 0.0
        vector<double> previousPz(config.numOfPoints + 1, 0.0); // Initialize to 0.0

        auto multipolarfilename = generate_multconvan_path(config);
        std::ofstream mout(multipolarfilename);
        mout.precision(17);
        mout << "Lmax\t\t\tAverage Relative Error\n";

        //Initialize Scattering Functions AA and BB when vv once vv stops changing for LSmax
        Initialize_ScatteringFunctions(AA, BB, config.velocity);

        // Start the timer
        auto start = std::chrono::high_resolution_clock::now();

        // ********************************************************************************
        // Start loop for Lmax
        // ********************************************************************************
        for (config.Lmax = config.minLmax; config.Lmax <= LSmax; config.Lmax++){

            // Print dynamic line
            cout << "Progress: " 
                  << "Lmax : " << config.Lmax
                  << " | v : " << config.velocity << "*c"
                  << " | NP radius: " << config.a << " nm" << endl << endl;

            auto filenamex = generate_LMT_path(config, "x");
            std::ofstream outx(filenamex);
            outx.precision(17);
            outx << "vv\tb\tDPExInt\t\t\tDPHxInt\t\t\tDPExScat\t\t\tDPHxScat\t\t\tDPExExt\t\t\tDPHxExt\t\t\tDPx\n";
            
            auto errorfilenamex = generate_LMT_path(config, "x", "error_");
            std::ofstream eoutx(errorfilenamex);
            eoutx.precision(17);
            eoutx << "vv\tb\terrDPExInt\t\t\terrDPHxInt\t\t\terrDPExScat\t\t\terrDPHxScat\t\t\terrDPExExt\t\t\terrDPHxExt\n";

            auto filenamez = generate_LMT_path(config, "z");
            std::ofstream outz(filenamez);
            outz.precision(17);
            outz << "vv\tb\tDPEzInt\t\t\tDPHzInt\t\t\tDPEzScat\t\t\tDPHzScat\t\t\tDPEzExt\t\t\tDPHzExt\t\t\tDPz\n";
            
            auto errorfilenamez = generate_LMT_path(config, "z", "error_");
            std::ofstream eoutz(errorfilenamez);
            eoutz.precision(17);
            eoutz << "vv\tb\terrDPEzInt\t\t\terrDPHzInt\t\t\terrDPEzScat\t\t\terrDPHzScat\t\t\terrDPEzExt\t\t\terrDPHzExt\n";

            double totalRelativeError = 0.0; // Accumulator for relative errors
            bool thresholdReached = false;

            for (int i = 0; i <= config.numOfPoints; ++i)
            {   
                config.b = (config.bInit 
                    + (config.bFin-config.bInit)*i/config.numOfPoints);

                printf("b : %4.1f   -->   ", config.b);
                DP(config, DPx, DPz);
                double Px=0.0;
                double Pz=0.0;

                for (int rr = 0; rr < 6; ++rr){ 
                Px += DPx[rr].real();
                Pz += DPz[rr].real();}

                // Calculate the error
                double errorPx = fabs(Px - previousPx[i]);
                double errorPz = fabs(Pz - previousPz[i]);
                double relativeError = (0.5*errorPx / fabs(Px))+(0.5*errorPz / fabs(Pz));

                // Accumulate the relative error
                totalRelativeError += relativeError;

               // Here print the total momentum
                outx << config.velocity << '\t'
                     << config.b << '\t'
                     << DPx[0].real() << '\t'
                     << DPx[1].real() << '\t'
                     << DPx[2].real() << '\t'
                     << DPx[3].real() << '\t'
                     << DPx[4].real() << '\t'
                     << DPx[5].real() << '\t'
                     << Px << '\n';
                eoutx << config.velocity << '\t'
                      << config.b << '\t'
                      << DPx[0].imag() << '\t'
                      << DPx[1].imag() << '\t'
                      << DPx[2].imag() << '\t'
                      << DPx[3].imag() << '\t'
                      << DPx[4].imag() << '\t'
                      << DPx[5].imag() << '\n';
                outz << config.velocity << '\t'
                     << config.b << '\t'
                     << DPz[0].real() << '\t'
                     << DPz[1].real() << '\t'
                     << DPz[2].real() << '\t'
                     << DPz[3].real() << '\t'
                     << DPz[4].real() << '\t'
                     << DPz[5].real() << '\t'
                     << Pz << '\n';
                eoutz << config.velocity << '\t'
                      << config.b << '\t'
                      << DPz[0].imag() << '\t'
                      << DPz[1].imag() << '\t'
                      << DPz[2].imag() << '\t'
                      << DPz[3].imag() << '\t'
                      << DPz[4].imag() << '\t'
                      << DPz[5].imag() << '\n';

                // Update previousP for the next iteration
                previousPx[i] = Px;
                previousPz[i] = Pz;
            }

            // Calculate the average relative error
            double averageRelativeError = totalRelativeError / (config.numOfPoints + 1);
            mout << config.Lmax << '\t' << averageRelativeError << '\n';
            cout << "Average Relative Error is: " << averageRelativeError << endl;
            cout << "------------------------------------------------------" << endl << endl;

            // Check if the average error is below the threshold
            if (averageRelativeError < config.errorThreshold) {
                thresholdReached = true;
                mout << "Convergence achieved at Lmax = " << config.Lmax
                     << " with average relative error = " << averageRelativeError << '\n';
                cout << "Convergence achieved at Lmax = " << config.Lmax 
                     << " with average relative error = " << averageRelativeError << endl;
                cout << "------------------------------------------------------" << endl << endl;
            }

            // Close the files
            outx.close();
            eoutx.close();
            outz.close();
            eoutz.close();

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
        cout << "Execution time: " << duration.count() << " seconds." << std::endl;
        mout << "Execution time: " << duration.count() << " seconds." << std::endl;
        mout.close();

        /*cout << "a = " << a / nm << "nm." << endl;
        cout << "r = " << r / nm << "nm." << endl;
        cout << "v = " << vv << "c." << endl;
        cout << "bInit = " << bInit / nm << "nm." << endl;
        cout << "bFin = " << bFin / nm << "nm." << endl;
        cout << endl;*/
    }
    else if (config.isBvsVContour){
        true;
    }
    else{

        // Store P values for each i in the previous Lmax iteration
        double previousPx = 0.0; // Initialize to 0.0
        double previousPz = 0.0; // Initialize to 0.0

        auto multipolarfilename = generate_multconvan_path(config);
        std::ofstream mout(multipolarfilename);
        mout.precision(17);
        mout << "Lmax\t\t\tAverage Relative Error\n";

        //Initialize Scattering Functions AA and BB when vv once vv stops changing for LSmax
        Initialize_ScatteringFunctions(AA, BB, config.velocity);

        // Start the timer
        auto start = std::chrono::high_resolution_clock::now();

        // ********************************************************************************
        // Start loop for Lmax
        // ********************************************************************************
        for (config.Lmax = config.minLmax; config.Lmax <= LSmax; config.Lmax++){

            // Print dynamic line
            cout << "Progress: " 
                  << "Lmax : " << config.Lmax
                  << " | v : " << config.velocity << "*c"
                  << " | b : " << config.b << " nm"
                  << " | NP radius: " << config.a << " nm" << endl << endl;

            auto filenamex = generate_LMT_path(config, "x");
            std::ofstream outx(filenamex);
            outx.precision(17);
            outx << "vv\tb\tDPExInt\t\t\tDPHxInt\t\t\tDPExScat\t\t\tDPHxScat\t\t\tDPExExt\t\t\tDPHxExt\t\t\tDPx\n";
            
            auto errorfilenamex = generate_LMT_path(config, "x", "error_");
            std::ofstream eoutx(errorfilenamex);
            eoutx.precision(17);
            eoutx << "vv\tb\terrDPExInt\t\t\terrDPHxInt\t\t\terrDPExScat\t\t\terrDPHxScat\t\t\terrDPExExt\t\t\terrDPHxExt\n";

            auto filenamez = generate_LMT_path(config, "z");
            std::ofstream outz(filenamez);
            outz.precision(17);
            outz << "vv\tb\tDPEzInt\t\t\tDPHzInt\t\t\tDPEzScat\t\t\tDPHzScat\t\t\tDPEzExt\t\t\tDPHzExt\t\t\tDPz\n";
            
            auto errorfilenamez = generate_LMT_path(config, "z", "error_");
            std::ofstream eoutz(errorfilenamez);
            eoutz.precision(17);
            eoutz << "vv\tb\terrDPEzInt\t\t\terrDPHzInt\t\t\terrDPEzScat\t\t\terrDPHzScat\t\t\terrDPEzExt\t\t\terrDPHzExt\n";

            bool thresholdReached = false;

            DP(config, DPx, DPz);
            double Px=0.0;
            double Pz=0.0;

            for (int rr = 0; rr < 6; ++rr){ 
                Px += DPx[rr].real();
                Pz += DPz[rr].real();}

            // Calculate the error
            double errorPx = fabs(Px - previousPx);
            double errorPz = fabs(Pz - previousPz);
            double relativeError = (0.5*errorPx / fabs(Px))+(0.5*errorPz / fabs(Pz));

            // Here print the total momentum
            outx << config.velocity << '\t'
                 << config.b << '\t'
                 << DPx[0].real() << '\t'
                 << DPx[1].real() << '\t'
                 << DPx[2].real() << '\t'
                 << DPx[3].real() << '\t'
                 << DPx[4].real() << '\t'
                 << DPx[5].real() << '\t'
                 << Px << '\n';
            eoutx << config.velocity << '\t'
                  << config.b << '\t'
                  << DPx[0].imag() << '\t'
                  << DPx[1].imag() << '\t'
                  << DPx[2].imag() << '\t'
                  << DPx[3].imag() << '\t'
                  << DPx[4].imag() << '\t'
                  << DPx[5].imag() << '\n';
            outz << config.velocity << '\t'
                 << config.b << '\t'
                 << DPz[0].real() << '\t'
                 << DPz[1].real() << '\t'
                 << DPz[2].real() << '\t'
                 << DPz[3].real() << '\t'
                 << DPz[4].real() << '\t'
                 << DPz[5].real() << '\t'
                 << Pz << '\n';
            eoutz << config.velocity << '\t'
                  << config.b << '\t'
                  << DPz[0].imag() << '\t'
                  << DPz[1].imag() << '\t'
                  << DPz[2].imag() << '\t'
                  << DPz[3].imag() << '\t'
                  << DPz[4].imag() << '\t'
                  << DPz[5].imag() << '\n';

            // Update previousP for the next iteration
            previousPx = Px;
            previousPz = Pz;

            // Calculate the average relative error
            mout << config.Lmax << '\t' << relativeError << '\n';
            cout << "Average Relative Error is: " << relativeError << endl;
            cout << "------------------------------------------------------" << endl << endl;

            // Check if the average error is below the threshold
            if (relativeError < config.errorThreshold) {
                thresholdReached = true;
                mout << "Convergence achieved at Lmax = " << config.Lmax
                     << " with average relative error = " << relativeError << '\n';
                cout << "Convergence achieved at Lmax = " << config.Lmax 
                     << " with average relative error = " << relativeError << endl;
                cout << "------------------------------------------------------" << endl << endl;
            }

            // Close the files
            outx.close();
            eoutx.close();
            outz.close();
            eoutz.close();

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
        cout << "Execution time: " << duration.count() << " seconds." << std::endl;
        mout << "Execution time: " << duration.count() << " seconds." << std::endl;
        mout.close();
    }

        
} // end void
