#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include "hf_aux.h"

using namespace std;

typedef std::vector<vector<double>> Real_Matrix;

int main(int argc, char* argv[])
{
//Assuming restricted hartree fock, read in arguments from commandline

   int ao = atoi(argv[1]);
   cout <<"Number of orbitals/el: " << ao << endl;

//Initializing energy
   double hf_energy=0;

//read in nuclear energy term
   double nuclear_en;
   read_nuc_en(nuclear_en);
    std::cout << "Nuclear energy is: " << nuclear_en << std::endl;

//read in kinetic energy terms
   cout << "Reading in Kinetic Energy Matrix" << endl;
   Real_Matrix T_int(ao-1,vector<double>(ao-1,0.0));
   read_T(ao, T_int);
   //double T_int[ao-1][ao-1];
   //T_int[ao-1][ao-1]=read_T(ao);
   
//read in overlap matrix
   cout << "Reading in Overlap Matrix" << endl;
   Real_Matrix S(ao-1,vector<double>(ao-1,0.0));
   read_S(ao,S);
   //double S[ao-1][ao-1];
   //S[ao-1][ao-1]=read_S(ao);
 
//read in v_int
   cout << "Reading in Interaction Matrix Matrix" << endl;
   Real_Matrix v_int(ao-1,vector<double>(ao-1,0.0));
   read_v_int(ao,v_int);
   //double v_int[ao-1][ao-1];
   //v_int[ao-1][ao-1]=read_v_int(ao);


   cout << "The energy is: " << hf_energy << endl;

   return 0;
}
