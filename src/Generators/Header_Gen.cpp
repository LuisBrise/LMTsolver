/*
 * Program: Header_Gen
 * Description:
 *   This program generates a header file (e.g., IN31.h) by combining a predefined template
 *   with the contents of two input files (e.g., FUNCTIONS_I31_2.dat and FUNCTIONS_N31_2.dat).
 *   The program dynamically names the output and input files based on a user-provided number N.
 *   For example, if N = 31, the program will:
 *     1. Create a file named IN31.h.
 *     2. Use FUNCTIONS_I31_2.dat and FUNCTIONS_N31_2.dat as input files.
 *     3. Write the combined content to IN31.h, including the template, the contents of the
 *        input files, and a footer.
 *
 *   The program is useful for automating the creation of header files with dynamic content
 *   based on a parameter N.
 *
 * How to Run:
 *   1. Compile the program using a C++ compiler (e.g., g++):
 *      g++ -o Header_Gen.out Header_Gen.cpp
 *   2. Run the program with the desired value of N as a command-line argument:
 *      ./Header_Gen.out <N>
 *      Example:
 *      ./Header_Gen.out 31
 *   3. The program will generate a file named IN<N>.h (e.g., IN31.h) in the same directory.
 *   4. Ensure that the input files (FUNCTIONS_I<N>_2.dat and FUNCTIONS_N<N>_2.dat) exist
 *      in the same directory as the program.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

int main(int argc, char* argv[]) {
    // Check if the user provided the number N as a command-line argument
    /*argc is an integer that represents the number of command-line arguments passed to the program.
    *It always has a value of at least 1, 
    *because the first argument (argv[0]) is always the name of the program itself.
    *
    *argv is an array of C-style strings (char*), where each element is a command-line argument.
    *argv[0] is the name of the program.
    *argv[1] to argv[argc-1] are the additional arguments provided by the user.
    */
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <N>" << std::endl;
        return 1;
    }

    // Parse the number N from the command-line argument
    int N;
    try {
        N = std::stoi(argv[1]);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: Invalid number N provided." << std::endl;
        return 1;
    }

    // Generate filenames based on N
    std::stringstream ss;
    ss << "IN" << N << ".h";
    std::string outputFileName = ss.str();

    ss.str(""); // Clear the stringstream
    ss << "FUNCTIONS_I" << N << "_2.dat";
    std::string inputFile1Name = ss.str();

    ss.str(""); // Clear the stringstream
    ss << "FUNCTIONS_L" << N << "_2.dat";
    std::string inputFile2Name = ss.str();

    // Open the output file
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not create " << outputFileName << " file!" << std::endl;
        return 1;
    }

    // Write the header content
    outFile << R"(#ifndef IN)" << N << R"(_H
#define IN)" << N << R"(_H

#include <iostream>                                                // Standart i/o C++ library
#include <fstream>
#include <cmath>
#include <chrono>
#include <complex>                                                 // Compĺex numbers
#include <boost/math/special_functions/bessel.hpp>                 // BOOST LIBRARIES:  1. BesselK in external fields
#include <boost/math/special_functions/beta.hpp>                   // Beta function for recursive relation in scatterred fields
#include <boost/math/special_functions/legendre.hpp>               // Lengendre Plm
#include <boost/math/quadrature/gauss_kronrod.hpp> 

// Spherical Bessel (Hankel) funcions with *complex argument* (very important) for the scatterred fields. 
// Only found in SLATEC project (FORTRAN 77).

#include <complex_bessel_bits/sph_besselFunctions.h>  

#define LSmax )" << N << R"(         //Maximun pole Lmax written in this file!

// Define the Integrals, using the notation IN(l1,l2,m1,m2) -> IN[l1-1][l2-1][m1+l1][m2+l2]

void FUN(
    double IM[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1],
    double IN[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1],
    double IU[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1],
    double IV[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1],
    double IW[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1],
    double IX[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1],
    double IY[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1],
    double IZ[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1],
    double ID[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1],
    double III[LSmax][LSmax+1][LSmax+1]){
)";

    // Write 10 line breaks
    for (int i = 0; i < 10; ++i) {
        outFile << std::endl;
    }

    // Open and write the content of FUNCTIONS_I<N>_2.dat
    std::ifstream inFile1(inputFile1Name);
    if (!inFile1.is_open()) {
        std::cerr << "Error: Could not open " << inputFile1Name << "!" << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(inFile1, line)) {
        outFile << line << std::endl;
    }
    inFile1.close();

    // Write 10 line breaks
    for (int i = 0; i < 10; ++i) {
        outFile << std::endl;
    }

    // Open and write the content of FUNCTIONS_N<N>_2.dat
    std::ifstream inFile2(inputFile2Name);
    if (!inFile2.is_open()) {
        std::cerr << "Error: Could not open " << inputFile2Name << "!" << std::endl;
        return 1;
    }

    while (std::getline(inFile2, line)) {
        outFile << line << std::endl;
    }
    inFile2.close();

    // Write the footer content
    outFile << R"(
}
#endif
)";

    // Close the output file
    outFile.close();

    std::cout << "File " << outputFileName << " has been created successfully." << std::endl;
    return 0;
}