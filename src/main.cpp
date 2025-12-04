// g++ -Iinclude -Iinclude/Dielectric/DrudeAl -c src/main.cpp -o build/main.o -lcomplex_bessel -w && g++ build/IN51.o build/main.o -o DrudeAl_LMTsolver_L51.out -lcomplex_bessel -w

#include "LMT.h"
//*********************************************************************************************************
//********************************************** MAIN PROGRAM *********************************************
int main(int argc, char** argv){
    // Default parameters
    double a, b;    // in nm
    double v;     // velocity (c units)
    bool isVScan = false; //Boolean for scan
    bool isBScan = false;
    bool isBvsVContour = false; // Boolean for contour 
    double threshold;
    int numOfPoints, minLmax;

    double bInit = std::numeric_limits<double>::quiet_NaN();
    double bFin  = std::numeric_limits<double>::quiet_NaN();
    double vInit = std::numeric_limits<double>::quiet_NaN();
    double vFin  = std::numeric_limits<double>::quiet_NaN();

    CLI::App app{"Linear Momentum Transfer"};

    app.add_option("-a", a, "Radius of the nanoparticle in nanometers [nm]")->default_val(1.0);
    app.add_option("-b", b, "Impact parameter (distance from the beam to the center of nanoparticle) in nanometers [nm]")->default_val(1.5);
    app.add_option("-v", v, "Electron velocity (as a fraction of the speed of light, c)")->default_val(0.7)->check(CLI::Range(0.0, 0.999));
    
    app.add_flag("--vscan", isVScan, "Enable scan over different electron velocities");
    app.add_option("--vinit", vInit, "Initial value of speed in scan");
    app.add_option("--vfin", vFin, "Final value of speed in scan");
    
    app.add_flag("--bscan", isBScan, "Enable scan over different impact parameters");
    app.add_option("--binit", bInit, "Initial value of impact parameter in scan");
    app.add_option("--bfin", bFin, "Final value of impact parameter in scan");
    
    app.add_flag("--contour", isBvsVContour, "Enable contour plot scan of angular momentum transfer as a function of velocity and impact parameter");
    
    app.add_option("-e", threshold, "Convergence threshold for the multipolar convergence")->default_val(pow(10, -4));
    app.add_option("-n", numOfPoints, "Number of sampling points for the scan or contour")->default_val(9);
    app.add_option("-m", minLmax, "Minimum multipole order (ℓ) for starting the multipolar convergence analysis")->default_val(1);

    CLI11_PARSE(app, argc, argv);

    // Validación
    if (b <= a) {
        cout << "Error: Impact parameter b (" << b << ") must be greater than NP radius (" << a << ")." << std::endl;
        cout << "Assigning value of " << a + 1.0 << " nm to b" << endl;
        b = a + 1.0;
    }

    // Assign default values after parsing
    if (isBScan) {
        if (std::isnan(bInit)) bInit = a + 0.5;
        if (std::isnan(bFin))  bFin  = a + 10.5;
    }
    if (isVScan) {
        if (std::isnan(vInit)) vInit = 0.50;
        if (std::isnan(vFin))  vFin  = 0.95;
    }
    if (isBvsVContour) {
        if (std::isnan(bInit)) bInit = a + 0.5;
        if (std::isnan(bFin))  bFin  = a + 10.5;
        if (std::isnan(vInit)) vInit = 0.50;
        if (std::isnan(vFin))  vFin  = 0.95;
    }

    cout.precision(17);    

    // 1-Define Physical parameters
/*  config(material, LSmax, a, b,vv,
           isVScan, isBScan, isBvsVContour,
           errorThreshold, numOfPoints, minLmax);*/
    SimulationConfig config(material,      
                            LSmax,          //Lmax
                            double(a),      //Np radius [nm]
                            double(b),      //Impact parameter [nm]
                            double(v),      // velocity [c units]
                            isVScan,
                            isBScan,
                            isBvsVContour,
                            bInit,
                            bFin,
                            vInit,
                            vFin,
                            threshold,     // error threshold
                            numOfPoints,              // number of points in the scan or contour
                            minLmax,              // min value of Lmax for multipolar convergence
                            get_current_timestamp()); 
    // 2-Calls the irreducibe integrals
    Initialize_Integrals(integrals);

    // 3-Calls the Gauss - Konrod function
    Omegas(xi, xk, xg);

    LMTsolver(config);
    //dcomplex DPx[6]; 
    //dcomplex DPz[6];
    //DP(config, DPx, DPz);    
 
return(0);

}
