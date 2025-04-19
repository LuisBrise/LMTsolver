/*
 * Program: II_Gen_Normalized
 * Description:
 *   This program calculates integrals involving Legendre polynomials and writes the results to a file.
 *   The integrals are computed using the Gauss-Kronrod quadrature method from the Boost library.
 *   The program dynamically generates a filename based on a user-provided value of Lmax, which determines
 *   the maximum order of the Legendre polynomials and the range of indices for the integrals.
 *
 *   The results are written to a file named FUNCTIONS_I<Lmax>_2.dat, where <Lmax> is the value provided by the user.
 *   Each line in the output file contains the result of an integral in the format:
 *     III[l-1][m][j] = <value>;
 *
 *   The program is useful for numerical analysis and scientific computing tasks involving Legendre polynomials.
 *
 * How to Run:
 *   1. Compile the program using a C++ compiler (e.g., g++):
 *      g++ -o II_Gen_Normalized.out II_Gen_Normalized.cpp -lboost_math_tr1
 *   2. Run the program:
 *      ./II_Gen_Normalized.out
 *   3. Enter the value of Lmax when prompted:
 *      Enter the value of Lmax: <your_value>
 *   4. The program will generate a file named FUNCTIONS_I<Lmax>_2.dat in the same directory.
 *   5. Ensure that the Boost library is installed and linked correctly during compilation.
 */

#include <stdio.h>
#include <boost/math/special_functions/beta.hpp>                   // Beta function for recursive relation in scatterred fields
#include <boost/math/special_functions/legendre.hpp>          
#include <boost/math/quadrature/gauss_kronrod.hpp>   
#include <fstream>          

double Pi     = boost::math::constants::pi<double>();

double LP(int l, int m, double x){return boost::math::legendre_p(l,m,x);}
//double L2P(int l1, int l2, int m1, int m2, double x){return LP(l1,m1,x)*LP(l2,m2,x);}

double factorial(int n){
    if (n == 0 || n==1){
      return 1.;
    }else{
    return n*factorial(n-1);
    }
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

double alm(int l, int m){
	if (abs(m)<=l)
	{
		return pow((2.*l+1.)*factorial(l-m)/(4.*Pi*factorial(l+m)),0.5);
	}
	else{return 0.;}
}

double II(int l, int m, int i1, int i2){

if (l == m - 2 || l == m - 1){return 0.;} 
     else if(l == m){ 	if (i2 % 2 == 0){
     		               return pow(-1.,m)*factorial2(2*m - 1)*boost::math::beta((i1 + m + 2)/2., (i2 + 1)/2.,1.);}
     		            else{return 0.;}}
     	else{return (1./(l - m))*((2.*l - 1.)*II(l - 1, m,i1, i2 + 1) - (l + m - 1.)*II(l - 2, m,i1, i2));}
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

char filename[sizeof "FUNCTIONS_I30_2.dat"];
sprintf(filename, "FUNCTIONS_I%d_2.dat",Lmax);
FILE *fpp = fopen(filename,"w+");

 if (!fpp) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return 1;
    }

double III;
for (int l = 1; l <= Lmax; ++l){
    for (int m = 0; m <= l; ++m){
        for (int j = m; j <= l; ++j){

        // We do not use the Condon Shortley Phase because in 1999 Garcia de Abajo's definition
        // of LegendreP did not include the phase. That is why we put pow(1.0,m) to remember that fact
        auto FI = [l,m,j](double x) { return alm(l, m)*pow(1.0,m)*pow(1-x*x,j/2.0)*pow(x,l-j)*LP(l,m,x);};

        III= boost::math::quadrature::gauss_kronrod<double, 181>::integrate( FI , -1., 1., 10, 1e-15, &error);
        fprintf(fpp,"III[%d][%d][%d] = %.18g;\n", l - 1, m, j, III);
            
        }
    }
    printf(" l = %d \n", l);
}
fclose(fpp);
std::cout << "Results written to " << filename << std::endl;

return 0;
}