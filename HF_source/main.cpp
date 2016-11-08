#include <iostream>
#include <fstream>
#include <iomanip>
//#include <cstdio>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "hf_aux.cpp"

using namespace std;

typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;


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
   Real_Matrix T_int(ao,vector<double>(ao,0.0));
   read_T(ao, T_int);


//read in overlap matrix
   cout << "Reading in Overlap Matrix" << endl;
   Real_Matrix S(ao,vector<double>(ao,0.0));
   read_S(ao,S);
 
//read in v_nuc
   cout << "Reading in Interaction Matrix Matrix" << endl;
   Real_Matrix v_nuc(ao,vector<double>(ao,0.0));
   read_v_nuc(ao,v_nuc);


//Build H_core = T_int + v_nuc
   cout << "Building H_core" << endl;
   Real_Matrix H_core(ao,vector<double>(ao,0.0));
   build_H_core(ao, v_nuc, T_int, H_core);

//Read v_int
   cout << "Reading in interaction matrix" << endl;
   Real_4dMatrix v_int(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
   read_v_int(ao, v_int);
   cout << "The energy is: " << hf_energy << endl;

   return 0;
}
