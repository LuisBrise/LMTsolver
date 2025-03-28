/*
 * Program: IN_Gen_Normalized
 * Description:
 *   This program calculates integrals involving Legendre polynomials and writes the results to a file.
 *   The integrals are computed using the Gauss-Kronrod quadrature method from the Boost library.
 *   The program dynamically generates a filename based on a user-provided value of Lmax, which determines
 *   the maximum order of the Legendre polynomials and the range of indices for the integrals.
 *
 *   The results are written to a file named FUNCTIONS_I<Lmax>_2.dat, where <Lmax> is the value provided by the user.
 *   Each line in the output file contains the result of an integral in the format:
 *     IM[l1 - 1][l2 -1][m1 + l1][m2 + l2] = <value>; etc...
 *
 *   The program is useful for numerical analysis and scientific computing tasks involving Legendre polynomials.
 *
 * How to Run:
 *   1. Compile the program using a C++ compiler (e.g., g++):
 *      g++ -o IN_Gen_Normalized.out IN_Gen_Normalized.cpp -lboost_math_tr1
 *   2. Run the program:
 *      ./IN_Gen_Normalized.out
 *   3. Enter the value of Lmax when prompted:
 *      Enter the value of Lmax: <your_value>
 *   4. The program will generate a file named FUNCTIONS_I<Lmax>_2.dat in the same directory.
 *   5. Ensure that the Boost library is installed and linked correctly during compilation.
 */

#include <stdio.h>
#include <boost/math/special_functions/legendre.hpp>          
#include <boost/math/quadrature/gauss_kronrod.hpp>   
#include <fstream>             

double Pi     = boost::math::constants::pi<double>();

double LP(int l, int m, double x){return boost::math::legendre_p(l,m,x);}
double L2P(int l1, int l2, int m1, int m2, double x){return LP(l1,m1,x)*LP(l2,m2,x);}

double factorial(int n){
    if (n == 0 || n==1){
      return 1.;
    }else{
    return n*factorial(n-1);
    }
}

double alm(int l, int m){
	if (abs(m)<=l)
	{
		return pow((2.*l+1.)*factorial(l-m)/(4.*Pi*factorial(l+m)),0.5);
	}
	else{return 0.;}
}

double IntD(int l1, int l2, int m1, int m2){
    if (l1 == l2 && m1 == m2){
        return 2.*alm(l1,m1)*alm(l2,m2)*factorial(l1+m1)/((2.*l1 + 1.)*factorial(l1-m1));}
        else{return 0.;}
}

int main(void){

int Lmax;
    double error;

    // Prompt the user to input Lmax
    std::cout << "Enter the value of Lmax: ";
    std::cin >> Lmax;

    // Validate the input
    if (Lmax <= 0) {
        std::cerr << "Error: Lmax must be a positive integer." << std::endl;
        return 1;
    }

char filename[sizeof "FUNCTIONS_L30_2.dat"];
sprintf(filename, "FUNCTIONS_L%d_2.dat",Lmax);
FILE *fpp = fopen(filename,"w+");

 if (!fpp) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return 1;
    }

double IM, IN, IU, IV, IW, IX, IY, IZ, ID;
for (int l1 = 1; l1 <= Lmax; ++l1){
    for (int l2 = 1; l2 <= Lmax; ++l2){
        for (int m1 = -l1; m1 <= l1; ++m1){
            for (int m2 = -l2; m2 <= l2; ++m2){
                auto FN = [l1,l2,m1,m2](double x) { return alm(l1,m1)*alm(l2,m2)*L2P(l1,l2,m1,m2,x)*sqrt(1.-x*x);};
                auto FM = [l1,l2,m1,m2](double x) { return alm(l1,m1)*alm(l2,m2)*L2P(l1,l2,m1,m2,x)*x;};
                auto FU = [l1,l2,m1,m2](double x) { return alm(l1,m1)*alm(l2,m2)*L2P(l1,l2,m1,m2,x)/sqrt(1.-x*x);};
                auto FV = [l1,l2,m1,m2](double x) { return alm(l1,m1)*alm(l2,m2)*L2P(l1,l2,m1,m2,x)*x*x/sqrt(1.-x*x);};
                auto FW = [l1,l2,m1,m2](double x) { return alm(l1,m1)*alm(l2,m2)*L2P(l1,l2,m1,m2,x)*x/sqrt(1.-x*x);};
                auto FX = [l1,l2,m1,m2](double x) { return alm(l1,m1)*alm(l2,m2)*L2P(l1,l2,m1,m2,x)*x*x/(1.-x*x);};
                auto FY = [l1,l2,m1,m2](double x) { return alm(l1,m1)*alm(l2,m2)*L2P(l1,l2,m1,m2,x)*x/(1.-x*x);};
                auto FZ = [l1,l2,m1,m2](double x) { return alm(l1,m1)*alm(l2,m2)*L2P(l1,l2,m1,m2,x)*x*x*x/(1.-x*x);};

                if(m2 == m1+1 || m2 == m1-1){
                IN = boost::math::quadrature::gauss_kronrod<double, 181>::integrate( FN , -1., 1., 10, 1e-15, &error);
                IU = boost::math::quadrature::gauss_kronrod<double, 181>::integrate( FU , -1., 1., 10, 1e-15, &error);
                IV = boost::math::quadrature::gauss_kronrod<double, 181>::integrate( FV , -1., 1., 10, 1e-15, &error);
                IW = boost::math::quadrature::gauss_kronrod<double, 181>::integrate( FW , -1., 1., 10, 1e-15, &error);

                fprintf(fpp,"IN[%d][%d][%d][%d] = %.18g;\n", l1 - 1, l2 -1, m1 + l1 , m2 + l2, IN);
                fprintf(fpp,"IU[%d][%d][%d][%d] = %.18g;\n", l1 - 1, l2 -1, m1 + l1 , m2 + l2, IU);
                fprintf(fpp,"IV[%d][%d][%d][%d] = %.18g;\n", l1 - 1, l2 -1, m1 + l1 , m2 + l2, IV);
                fprintf(fpp,"IW[%d][%d][%d][%d] = %.18g;\n", l1 - 1, l2 -1, m1 + l1 , m2 + l2, IW);
                }

                else if(m2 == m1){
                IM = boost::math::quadrature::gauss_kronrod<double, 181>::integrate( FM , -1., 1., 10, 1e-15, &error);
                IX = boost::math::quadrature::gauss_kronrod<double, 181>::integrate( FX , -1., 1., 10, 1e-15, &error);
                IY = boost::math::quadrature::gauss_kronrod<double, 181>::integrate( FY , -1., 1., 10, 1e-15, &error);
                IZ = boost::math::quadrature::gauss_kronrod<double, 181>::integrate( FZ , -1., 1., 10, 1e-15, &error);
                ID = IntD(l1,l2,m1,m2);

                fprintf(fpp,"IM[%d][%d][%d][%d] = %.18g;\n", l1 - 1, l2 -1, m1 + l1 , m2 + l2, IM);
                fprintf(fpp,"IX[%d][%d][%d][%d] = %.18g;\n", l1 - 1, l2 -1, m1 + l1 , m2 + l2, IX);
                fprintf(fpp,"IY[%d][%d][%d][%d] = %.18g;\n", l1 - 1, l2 -1, m1 + l1 , m2 + l2, IY);
                fprintf(fpp,"IZ[%d][%d][%d][%d] = %.18g;\n", l1 - 1, l2 -1, m1 + l1 , m2 + l2, IZ);
                fprintf(fpp,"ID[%d][%d][%d][%d] = %.18g;\n", l1 - 1, l2 -1, m1 + l1 , m2 + l2, ID);
                }
            }
        }
    }
    printf(" l = %d \n", l1);
}
fclose(fpp);
std::cout << "Results written to " << filename << std::endl;

return 0;
}