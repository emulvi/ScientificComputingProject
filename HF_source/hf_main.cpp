#include <iostream>
#include <fstream>
#include <iomanip>
//#include <cstdio>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <string.h>
#include <stdio.h>
//#include <lapacke.h>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Cholesky"
#include "Eigen/Core"
#include "hf_aux.cpp"
//#include "mp2_aux.cpp"

using namespace std;
//using namespace Eigen;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;



double hf_main(int ao, int occ, char *Path){
//Assuming restricted hartree fock, read in arguments from commandline

//   int ao;
//   int occ;
//   char* Path;

   //if (i + 1 != argc){ // Check that we haven't finished parsing already

//Initializing energy
   double hf_energy=0;

//read in nuclear energy term
   double En_nuc = 0;
   En_nuc = read_nuc_en(Path);

//read in kinetic energy terms
   cout << "Reading in Kinetic Energy Matrix" << endl;
   Matrix T_int = Matrix::Zero(ao,ao);
   //Real_Matrix T_int(ao,vector<double>(ao,0.0));
   read_T(ao, T_int,Path);


//read in overlap matrix
   Matrix S = Matrix::Zero(ao,ao);
   read_S(ao,S,Path);
 
//read in v_nuc
   Matrix v_nuc  = Matrix::Zero(ao,ao);
   read_v_nuc(ao,v_nuc,Path);


//Build H_core = T_int + v_nuc
   Matrix H_core = Matrix::Zero(ao,ao);
   build_H_core(ao, v_nuc, T_int, H_core);

//Read v_int
   Matrix V = Matrix::Zero(ao*ao,ao*ao);
   read_v_int(ao, V,Path);

//calculate orthogonalization matrix, S^{-1/2}
   Matrix S12 = Matrix::Zero(ao,ao);
   Matrix Xmat = Matrix::Zero(ao,ao);
   calculate_S12(ao, S, S12, Xmat);

//Diagonalize Fock
   Matrix Fock = Matrix::Zero(ao,ao);
   Matrix C_ao = Matrix::Zero(ao,ao);
   Matrix evals = Matrix::Zero(ao,ao);
   Matrix evecs = Matrix::Zero(ao,ao);
   diagonalize_Fock(ao, H_core, Xmat, C_ao, evals, evecs);

//Build density matrix


//Build new Fock matrix, call it G
   Matrix P0 = Matrix::Zero(ao,ao);
   build_P(ao, occ, C_ao, P0);


//initial SCF electronic energy


   Matrix Fock2 = H_core;

   double En_elec =calculate_En_elec(ao, P0, H_core, Fock2);
   double En_total = En_elec + En_nuc;




////Build new Fock matrix, call it G
//   cout << "Building new Fock matrix, G" << endl;
//   Matrix Fock_new = Matrix::Zero(ao,ao);
//   build_new_Fock(ao, P0, v_int, H_core, Fock_new);
   Matrix Fock_new = Matrix::Zero(ao,ao);
   build_new_Fock(ao, P0, V, H_core, Fock_new);

//Build new Density Matrix
   Matrix P = Matrix::Zero(ao,ao);
   Matrix C_ao_new = Matrix::Zero(ao,ao);
   diagonalize_Fock(ao, Fock_new, Xmat, C_ao_new, evals, evecs);
   build_P(ao, occ, C_ao_new, P);

   double En_elec_new =calculate_En_elec(ao, P, H_core, Fock_new);
   Matrix P2 = Matrix::Zero(ao,ao);
   P2=P;

   int iteration = 2;
////////////Begin SCF Procedure
   double deltaE = 10;
   Matrix C_mo = Matrix::Zero(ao,ao);
   while (abs(deltaE) > 0.000000000001) {
         En_elec=En_elec_new;
         P0=P2;
	 //cout << "Density matrix: \n" << P0 << endl;

      //Build new Fock matrix, call it G
         Matrix Fock_new = Matrix::Zero(ao,ao);
         build_new_Fock(ao, P0, V, H_core, Fock_new);
      
      //Build new Density Matrix
         Matrix C_ao_new = Matrix::Zero(ao,ao);
	 Matrix P = Matrix::Zero(ao,ao);
         diagonalize_Fock(ao, Fock_new, Xmat, C_ao_new, evals, evecs);
         build_P(ao, occ, C_ao_new, P);

      //Compute new SCF Energy
         En_elec_new =calculate_En_elec(ao, P, H_core, Fock_new);
         En_total = En_elec_new + En_nuc;
	 
	 iteration = iteration + 1;
      
         deltaE = En_elec_new-En_elec;
	 P2=P;
         cout << "deltaE is: " << deltaE << endl;
	 C_mo = C_ao_new;

   }
   return En_total;
}
