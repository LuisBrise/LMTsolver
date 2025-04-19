// g++ -DLSmax=11 main.cpp IN11.cpp -o run_integrals
#ifndef IN11_H
#define IN11_H

#define LSmax 11         //Maximun pole Lmax written in this file!      

struct Irreducible_Integrals{
 double IM[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
 double IN[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
 double IU[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
 double IV[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
 double IW[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
 double IX[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
 double IY[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
 double IZ[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
 double ID[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
 double III[LSmax][LSmax+1][LSmax+1];
};

// Define the Integrals, using the notation IN(l1,l2,m1,m2) -> IN[l1-1][l2-1][m1+l1][m2+l2]

void Initialize_Integrals(Irreducible_Integrals& integrals);

#endif
