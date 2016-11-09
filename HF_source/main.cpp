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

// Test Comment by Sahil to check git push 

/* DSYEV prototype */
extern void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );


typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;


int main(int argc, char* argv[])
{
//Assuming restricted hartree fock, read in arguments from commandline

   int ao = atoi(argv[1]);
   int occ = atoi(argv[2]);
   cout <<"Number of orbitals/el: " << ao << endl;

//Initializing energy
   double hf_energy=0;

//read in nuclear energy term
   double En_nuc = read_nuc_en();
   std::cout << "Nuclear energy is: " << En_nuc << std::endl;

//read in kinetic energy terms
   cout << "Reading in Kinetic Energy Matrix" << endl;
   Matrix T_int = Matrix::Zero(ao,ao);
   //Real_Matrix T_int(ao,vector<double>(ao,0.0));
   read_T(ao, T_int);


//read in overlap matrix
   cout << "Reading in Overlap Matrix" << endl;
   Matrix S = Matrix::Zero(ao,ao);
   read_S(ao,S);
 
//read in v_nuc
   cout << "Reading in Interaction Matrix Matrix" << endl;
   Matrix v_nuc  = Matrix::Zero(ao,ao);
   read_v_nuc(ao,v_nuc);


//Build H_core = T_int + v_nuc
   cout << "Building H_core" << endl;
   Matrix H_core = Matrix::Zero(ao,ao);
   build_H_core(ao, v_nuc, T_int, H_core);

//Read v_int
   cout << "Reading in interaction matrix" << endl;
   Real_4dMatrix v_int(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
   read_v_int(ao, v_int);

//calculate orthogonalization matrix, S^{-1/2}
   cout << "Calculating S^{-1/2}" << endl;
   Matrix S12 = Matrix::Zero(ao,ao);
   Matrix Xmat = Matrix::Zero(ao,ao);
   calculate_S12(ao, S, S12, Xmat);

//Diagonalize Fock
   cout << "Diagonalizing Fock" << endl;
   Matrix Fock = Matrix::Zero(ao,ao);
   Matrix C_ao = Matrix::Zero(ao,ao);
   diagonalize_Fock(ao, H_core, Xmat, Fock, C_ao);

//Build density matrix
   cout << "Building initial guess for Density Matrix" << endl;
   Matrix P0 = Matrix::Zero(ao,ao);
   build_P(ao, occ, C_ao, P0);

//initial SCF electronic energy

   cout << "calculating the initial SCF energy: " << endl;

   double En_elec =calculate_En_elec(ao, P0, H_core, Fock);
   double En_total = En_elec + En_nuc;

   cout << "Total Energy is ...." << endl;
   cout << En_elec << "+" << En_nuc << "=" << En_total << endl;


//Build new Fock matrix, call it G
   cout << "Building new Fock matrix, G" << endl;
   Matrix Fock_new = Matrix::Zero(ao,ao);
   build_new_Fock(ao, P0, v_int, H_core, Fock_new);

//Build new Density Matrix
   cout << "Building new Density matrix, Pnew" << endl;
   Matrix P = Matrix::Zero(ao,ao);
   diagonalize_Fock(ao, H_core, Xmat, Fock_new, C_ao);
   build_P(ao, occ, C_ao, P);

//Compute new SCF Energy
   En_elec =calculate_En_elec(ao, P, H_core, Fock_new);
   En_total = En_elec + En_nuc;

   cout << "Total Energy is ...." << endl;
   cout << En_elec << "+" << En_nuc << "=" << En_total << endl;


   cout << "The energy is: " << hf_energy << endl;

   return 0;
}



