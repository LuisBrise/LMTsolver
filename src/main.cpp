// g++ -Iinclude -Iinclude/Dielectric/DrudeAl -c src/main.cpp -o build/main.o -lcomplex_bessel -w && g++ build/IN51.o build/main.o -o DrudeAl_AMTsolver_L51.out -lcomplex_bessel -w

#include "AMT.h"
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

    CLI::App app{"Angular Momentum Transfer"};

    app.add_option("-a", a, "Radius of the nanoparticle in nanometers [nm]")->default_val(1.0);
    app.add_option("-b", b, "Impact parameter (distance from the beam to the center of nanoparticle) in nanometers [nm]")->default_val(1.5);
    app.add_option("-v", v, "Electron velocity (as a fraction of the speed of light, c)")->default_val(0.7)->check(CLI::Range(0.0, 0.999));
    
    app.add_flag("--vscan", isVScan, "Enable scan over different electron velocities");
    app.add_flag("--bscan", isBScan, "Enable scan over different impact parameters");
    app.add_flag("--contour", isBvsVContour, "Enable contour plot scan of angular momentum transfer as a function of velocity and impact parameter");
    
    app.add_option("-e", threshold, "Convergence threshold for the multipolar convergence")->default_val(pow(10, -4));
    app.add_option("-n", numOfPoints, "Number of sampling points for the scan or contour")->default_val(9);
    app.add_option("-m", minLmax, "Minimum multipole order (ℓ) for starting the multipolar convergence analysis")->default_val(1);

    CLI11_PARSE(app, argc, argv);

    // Validación
    if (b <= a) {
        std::cerr << "Error: Impact parameter b (" << b << ") must be greater than NP radius (" << a << ")." << std::endl;
        return 1;
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
                            threshold,     // error threshold
                            numOfPoints,              // number of points in the scan or contour
                            minLmax,              // min value of Lmax for multipolar convergence
                            get_current_timestamp()); 
    // 2-Calls the irreducibe integrals
    Initialize_Integrals(integrals);

    // 3-Calls the Gauss - Konrod function
    Omegas(xi, xk, xg);

    AMTsolver(config);
    //parallelAMTsolver(config);    
 
return(0);

}
