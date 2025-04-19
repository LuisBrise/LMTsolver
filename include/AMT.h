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
R"(AMTsolver

Description:
This program calculates the spectral angular momentum (dL/dω) transferred to a nanoparticle (NP) due to the interaction with a fast-moving electron. The calculation is performed by analytically solving the closed surface integral, resulting in a double multipolar sum expression. The program focuses exclusively on the external field contribution.This program calculates the total angular momentum transfer (DeltaL) for different values of the impact parameter b.

Key Features:
- The frequency integral is computed using the Gauss-Kronrod quadrature method, with an adaptive partition for rapid convergence.
- The spectral results are saved to files named "dldw_*.dat".
- The integrated angular momentum results are saved to files named "DL_*.dat".

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
     << "Minimum value of Lmax : " << config.minLmax << endl;
}

void DL(SimulationConfig& config, dcomplex DLy[6]){


double r = config.r*nm;
double vv = config.velocity;
double b = config.b*nm;
double a = config.a*nm;
int Lmax = config.Lmax;

double dm;
double dl;

double dldwy[6];

dcomplex CM[Lmax][2*Lmax +1], tMl[Lmax];
dcomplex DE[Lmax][2*Lmax +1], tEl[Lmax];
dcomplex dz[2][Lmax], df[2][Lmax];

dcomplex IErty[4], IHrty[4];
dcomplex IErfy[4], IHrfy[4];

double dl1, dl2, dm1, dm2;
double normfactor; // A normalization factor due to the normalization of integrals

double IN1, IV1, IW1, IW2, IU1, IU2;

double w, k0;

auto filename = generate_dldw_path(config, "y");
std::ofstream out(filename);
//write_metadata(filename.parent_path(), config);
out << "vv\t\t\tb\t\t\tw(au)\t\t\tdlEintdwy\t\t\tdlHintdwy\t\t\tdlEsdwy\t\t\tdlHsdwy\t\t\tdlEsdwyExt\t\t\tdlHsdwyExt\n";

dcomplex aa, bb;

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

                        IV1 = integrals.IV[l1-1][l2-1][m1+l1][m2+l2];

                        IW1 = integrals.IW[l1-1][l2-1][m1+l1][m2+l2];
                        IW2 = normfactor*integrals.IW[l1][l2-1][m1+l1+1][m2+l2];   

                        IU1 = integrals.IU[l1-1][l2-1][m1+l1][m2+l2];
                        IU2 = normfactor*integrals.IU[l1][l2-1][m1+l1+1][m2+l2];

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
out << vv << '\t'
         << b << '\t'
         << dldwy[0] << '\t'
         << dldwy[1] << '\t'
         << dldwy[2] << '\t'
         << dldwy[3] << '\t'
         << dldwy[4] << '\t'
         << dldwy[5] << '\n';

}// for w

out.close();

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

cout << "DLy : " << (DLy[0] + DLy[1] + DLy[2] + DLy[3]).real()
     << " +- " <<   (DLy[0] + DLy[1] + DLy[2] + DLy[3]).imag() << endl;
//cout << endl;

} //end void

void AMTsolver(SimulationConfig& config){

    print_initial_message(config);
    
    dcomplex DLy[6];

    std::string userResponse;
    bool validInput = false;

    if (config.isVScan){
        double vvInit = 0.50; 
        double vvFin  = 0.95;

        do{
            cout << "\nSuggested speed scan range:\n"
              << "  - Initial value (vInit): " << vvInit << "c\n"
              << "  - Final value (vFin):   " << vvFin  << " nm\n"
              << "Do you want to proceed with these values? [y/n]: ";

            std::getline(std::cin, userResponse);

            if (userResponse == "n" || userResponse == "no"){
                cout << "Select value for initial speed in range (0-0.99):";
                cin >> vvInit;
                cout << endl;
                cout << "Select value for final speed in range (0-0.99):";
                cin >> vvFin;
                cout << endl;
                validInput = true; 
            }
            else if (userResponse == "y" || userResponse == "yes"){
                validInput = true;
            }
            else {
                std::cout << "Invalid input. Please type 'y' or 'n'.\n";
            }
        } while (!validInput);

        // Vector to store Ly values for each i in the previous Lmax iteration
        vector<double> previousLy(config.numOfPoints + 1, 0.0); // Initialize to 0.0

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
            auto filename = generate_AMT_path(config, "y");
            std::ofstream out(filename);
            out.precision(17);
            out << "vv\tb\tDLEy\t\t\tDLHy\t\t\tDLEsy\t\t\tDLHsy\t\t\tDLEsyExt\t\t\tDLHsyExt\t\t\tDLy\n";
            
            auto errorfilename = generate_AMT_path(config, "y", "error_");
            std::ofstream eout(errorfilename);
            eout.precision(17);
            eout << "vv\tb\terrDLEy\t\t\terrDLHy\t\t\terrDLEsy\t\t\terrDLHsy\t\t\terrDLEsyExt\t\t\terrDLHsyExt\n";

            double totalRelativeError = 0.0; // Accumulator for relative errors
            bool thresholdReached = false;

            for (int i = 0; i <= config.numOfPoints; ++i)
            {
                config.velocity = (vvInit + (vvFin-vvInit)*i/config.numOfPoints);
                //Initialize Scattering Functions AA and BB when vv once vv stops changing for LSmax
                Initialize_ScatteringFunctions(AA, BB, config.velocity);
                printf("v : %3.2f   -->   ", config.velocity);
                DL(config, DLy);
                double Ly=0.0;

                for (int rr = 0; rr < 6; ++rr){ 
                    Ly += DLy[rr].real();}

                // Calculate the error
                double errorLy = fabs(Ly - previousLy[i]);
                double relativeError = errorLy / fabs(Ly);

                // Accumulate the relative error
                totalRelativeError += relativeError;

               // Here print the total momentum
               // Here print the total momentum
                out << config.velocity << '\t'
                    << config.b << '\t'
                    << DLy[0].real() << '\t'
                    << DLy[1].real() << '\t'
                    << DLy[2].real() << '\t'
                    << DLy[3].real() << '\t'
                    << DLy[4].real() << '\t'
                    << DLy[5].real() << '\t'
                    << Ly << '\n';
                eout << config.velocity << '\t'
                     << config.b << '\t'
                     << DLy[0].imag() << '\t'
                     << DLy[1].imag() << '\t'
                     << DLy[2].imag() << '\t'
                     << DLy[3].imag() << '\t'
                     << DLy[4].imag() << '\t'
                     << DLy[5].imag() << '\n';
                // Update previousLy for the next iteration
                previousLy[i] = Ly;
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
            out.close();
            eout.close();

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

        /*cout << "a = " << a/nm << "nm." << endl;
        cout << "r = " << r/nm << "nm." << endl;
        cout << "vvInit = " << vvInit << "c." << endl;
        cout << "vvFin = " << vvFin << "c." << endl;
        cout << "b = " << b/nm << "nm." << endl;
        cout << endl;*/
    }

    else if (config.isBScan)
    {
        double bInit = config.a + 0.5; 
        double bFin  = bInit + 10;

        do{
            cout << "\nSuggested impact parameter scan range:\n"
              << "  - Initial value (bInit): " << bInit << " nm\n"
              << "  - Final value (bFin):   " << bFin  << " nm\n"
              << "Do you want to proceed with these values? [y/n]: ";
            
            std::getline(std::cin, userResponse);

            if (userResponse == "n" || userResponse == "no"){
                cout << "Select value for initial impact parameter, greater than "<< config.a <<" :";
                cin >> bInit;
                cout << endl;
                cout << "Select value for final impact parameter, greater than "<< config.a <<" :";
                cin >> bFin;
                cout << endl;
                validInput = true; 
            } 
            else if (userResponse == "y" || userResponse == "yes"){
                validInput = true;
            }
            else{
                std::cout << "Invalid input. Please type 'y' or 'n'.\n";
            }
        } while (!validInput);

        // Vector to store Ly values for each i in the previous Lmax iteration
        vector<double> previousLy(config.numOfPoints + 1, 0.0); // Initialize to 0.0

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

            auto filename = generate_AMT_path(config, "y");
            std::ofstream out(filename);
            out.precision(17);
            out << "vv\tb\tDLEy\t\t\tDLHy\t\t\tDLEsy\t\t\tDLHsy\t\t\tDLEsyExt\t\t\tDLHsyExt\t\t\tDLy\n";
            
            auto errorfilename = generate_AMT_path(config, "y", "error_");
            std::ofstream eout(errorfilename);
            eout.precision(17);
            eout << "vv\tb\terrDLEy\t\t\terrDLHy\t\t\terrDLEsy\t\t\terrDLHsy\t\t\terrDLEsyExt\t\t\terrDLHsyExt\n";

            double totalRelativeError = 0.0; // Accumulator for relative errors
            bool thresholdReached = false;

            for (int i = 0; i <= config.numOfPoints; ++i)
            {   
                config.b = (bInit + (bFin-bInit)*i/config.numOfPoints);

                printf("b : %4.1f   -->   ", config.b);
                DL(config, DLy);
                double Ly=0.0;

                for (int rr = 0; rr < 6; ++rr){ 
                Ly += DLy[rr].real();}

                // Calculate the error
                double errorLy = fabs(Ly - previousLy[i]);
                double relativeError = errorLy / fabs(Ly);

                // Accumulate the relative error
                totalRelativeError += relativeError;

               // Here print the total momentum
                out << config.velocity << '\t'
                    << config.b << '\t'
                    << DLy[0].real() << '\t'
                    << DLy[1].real() << '\t'
                    << DLy[2].real() << '\t'
                    << DLy[3].real() << '\t'
                    << DLy[4].real() << '\t'
                    << DLy[5].real() << '\t'
                    << Ly << '\n';
                eout << config.velocity << '\t'
                     << config.b << '\t'
                     << DLy[0].imag() << '\t'
                     << DLy[1].imag() << '\t'
                     << DLy[2].imag() << '\t'
                     << DLy[3].imag() << '\t'
                     << DLy[4].imag() << '\t'
                     << DLy[5].imag() << '\n';

                // Update previousLy for the next iteration
                previousLy[i] = Ly;
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
            out.close();
            eout.close();

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

        // Store Ly values for each i in the previous Lmax iteration
        double previousLy = 0.0; // Initialize to 0.0

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

            auto filename = generate_AMT_path(config, "y");
            std::ofstream out(filename);
            out.precision(17);
            out << "vv\tb\tDLEy\t\t\tDLHy\t\t\tDLEsy\t\t\tDLHsy\t\t\tDLEsyExt\t\t\tDLHsyExt\t\t\tDLy\n";
            
            auto errorfilename = generate_AMT_path(config, "y", "error_");
            std::ofstream eout(errorfilename);
            eout.precision(17);
            eout << "vv\tb\terrDLEy\t\t\terrDLHy\t\t\terrDLEsy\t\t\terrDLHsy\t\t\terrDLEsyExt\t\t\terrDLHsyExt\n";

            bool thresholdReached = false;

            DL(config, DLy);
            double Ly=0.0;

            for (int rr = 0; rr < 6; ++rr){ 
                Ly += DLy[rr].real();}

            // Calculate the error
            double errorLy = fabs(Ly - previousLy);
            double relativeError = errorLy / fabs(Ly);

            // Here print the total momentum
            out << config.velocity << '\t'
                << config.b << '\t'
                << DLy[0].real() << '\t'
                << DLy[1].real() << '\t'
                << DLy[2].real() << '\t'
                << DLy[3].real() << '\t'
                << DLy[4].real() << '\t'
                << DLy[5].real() << '\t'
                << Ly << '\n';
            eout << config.velocity << '\t'
                 << config.b << '\t'
                 << DLy[0].imag() << '\t'
                 << DLy[1].imag() << '\t'
                 << DLy[2].imag() << '\t'
                 << DLy[3].imag() << '\t'
                 << DLy[4].imag() << '\t'
                 << DLy[5].imag() << '\n';

            // Update previousLy for the next iteration
            previousLy = Ly;

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
            out.close();
            eout.close();

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
