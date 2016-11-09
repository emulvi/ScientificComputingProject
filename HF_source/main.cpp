#include <iostream>
#include <fstream>
#include <iomanip>
//#include <cstdio>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
//#include <lapacke.h>
#include <mkl.h>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
#include "hf_aux.cpp"

using namespace std;

/* DSYEV prototype */
extern void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
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
   Matrix T_int(ao,ao);
   //Real_Matrix T_int(ao,vector<double>(ao,0.0));
   read_T(ao, T_int);


//read in overlap matrix
   cout << "Reading in Overlap Matrix" << endl;
   Matrix S(ao,ao);
   read_S(ao,S);
 
//read in v_nuc
   cout << "Reading in Interaction Matrix Matrix" << endl;
   Matrix v_nuc(ao,ao);
   read_v_nuc(ao,v_nuc);


//Build H_core = T_int + v_nuc
   cout << "Building H_core" << endl;
   Matrix H_core(ao,ao);
   build_H_core(ao, v_nuc, T_int, H_core);

//Read v_int
   cout << "Reading in interaction matrix" << endl;
   Real_4dMatrix v_int(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
   read_v_int(ao, v_int);

//calculate orthogonalization matrix, S^{-1/2}
   cout << "Calculating S^{-1/2}" << endl;
   Matrix S12(ao,ao);
   Matrix Xmat(ao,ao);
   calculate_S12(ao, S, S12, Xmat);

//Diagonalize Fock
   cout << "Diagonalizing Fock" << endl;
   Matrix Fock(ao,ao);
   Matrix C_ao(ao,ao);
   diagonalize_Fock(ao, H_core, Xmat, Fock, C_ao);

//Build density matrix
   cout << "Building initial guess for Density Matrix" << endl;
   Matrix P0(ao,ao);
   build_P(ao, C_ao, P0);

   cout << "The energy is: " << hf_energy << endl;

   return 0;
}
