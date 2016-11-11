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
#include "Eigen/Core"
#include "hf_aux.cpp"
#include "mp2_aux.cpp"

using namespace std;


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;


int main(int argc, char* argv[])
{
//Assuming restricted hartree fock, read in arguments from commandline

   int ao;
   int occ;
   char* Path;

   //if (i + 1 != argc){ // Check that we haven't finished parsing already
   for (int i=1; i<7; i=i+2){

     if (strncmp(argv[i],"-ao",3)==0) {
         ao =  atoi(argv[i + 1]);
     } 
     else if (strncmp(argv[i],"-occ",4)==0) {
       occ = atoi(argv[i + 1]);
     } 
     else if(strncmp(argv[i], "-path",5)==0) {
       Path = argv[i + 1];
     };
   };

   cout << "ao is: " << ao << endl;
   cout << "occ is " << occ << endl;
   cout << "Path is " << Path << endl;


   cout <<"Number of orbitals/el: " << ao << endl;

//Initializing energy
   double hf_energy=0;

//read in nuclear energy term
   double En_nuc = read_nuc_en(Path);
   std::cout << "Nuclear energy is: " << En_nuc << std::endl;

//read in kinetic energy terms
   cout << "Reading in Kinetic Energy Matrix" << endl;
   Matrix T_int = Matrix::Zero(ao,ao);
   //Real_Matrix T_int(ao,vector<double>(ao,0.0));
   read_T(ao, T_int,Path);


//read in overlap matrix
   cout << "Reading in Overlap Matrix" << endl;
   Matrix S = Matrix::Zero(ao,ao);
   read_S(ao,S,Path);
 
//read in v_nuc
   cout << "Reading in Interaction Matrix Matrix" << endl;
   Matrix v_nuc  = Matrix::Zero(ao,ao);
   read_v_nuc(ao,v_nuc,Path);


//Build H_core = T_int + v_nuc
   cout << "Building H_core" << endl;
   Matrix H_core = Matrix::Zero(ao,ao);
   build_H_core(ao, v_nuc, T_int, H_core);

//Read v_int
   cout << "Reading in interaction matrix" << endl;
   Real_4dMatrix v_int(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
   read_v_int(ao, v_int,Path);

//calculate orthogonalization matrix, S^{-1/2}
   cout << "Calculating S^{-1/2}" << endl;
   Matrix S12 = Matrix::Zero(ao,ao);
   Matrix Xmat = Matrix::Zero(ao,ao);
   calculate_S12(ao, S, S12, Xmat);

//Diagonalize Fock
   cout << "Diagonalizing Fock" << endl;
   Matrix Fock = Matrix::Zero(ao,ao);
   Matrix C_ao = Matrix::Zero(ao,ao);
   Matrix evals = Matrix::Zero(ao,ao);
   Matrix evecs = Matrix::Zero(ao,ao);
   diagonalize_Fock(ao, H_core, Xmat, C_ao, evals, evecs);

//Build density matrix
   cout << "Building initial guess for Density Matrix" << endl;


//Build new Fock matrix, call it G
   cout << "Building new Fock matrix, G" << endl;
   Matrix P0 = Matrix::Zero(ao,ao);
   build_P(ao, occ, C_ao, P0);


//initial SCF electronic energy

   cout << "calculating the initial SCF energy: " << endl;

   Matrix Fock2 = H_core;

   double En_elec =calculate_En_elec(ao, P0, H_core, Fock2);
   double En_total = En_elec + En_nuc;

   cout << "Total Energy is ...." << endl;
   cout << En_elec << "+" << En_nuc << "=" << En_total << endl;



////Build new Fock matrix, call it G
//   cout << "Building new Fock matrix, G" << endl;
//   Matrix Fock_new = Matrix::Zero(ao,ao);
//   build_new_Fock(ao, P0, v_int, H_core, Fock_new);
   Matrix Fock_new = Matrix::Zero(ao,ao);
   build_new_Fock(ao, P0, v_int, H_core, Fock_new);

//Build new Density Matrix
   cout << "Building new Density matrix, Pnew" << endl;
   Matrix P = Matrix::Zero(ao,ao);
   Matrix C_ao_new = Matrix::Zero(ao,ao);
   diagonalize_Fock(ao, Fock_new, Xmat, C_ao_new, evals, evecs);
   build_P(ao, occ, C_ao_new, P);

   double En_elec_new =calculate_En_elec(ao, P, H_core, Fock_new);
   Matrix P2 = Matrix::Zero(ao,ao);
   P2=P;

   cout << "before procedure En_elec new is : " << En_elec_new << endl;
   int iteration = 2;
////////////Begin SCF Procedure
   double deltaE = 10;
   Matrix C_mo = Matrix::Zero(ao,ao);
   while (abs(deltaE) > 0.000000000001) {
         cout << "--------------------------------------Iteration: " << iteration <<"-------------------------------------"<< endl;
         En_elec=En_elec_new;
	 cout << "En_elec is: " << En_elec << endl;
         P0=P2;
	 //cout << "Density matrix: \n" << P0 << endl;

      //Build new Fock matrix, call it G
         cout << "Building new Fock matrix, G" << endl;
         Matrix Fock_new = Matrix::Zero(ao,ao);
         build_new_Fock(ao, P0, v_int, H_core, Fock_new);
      
      //Build new Density Matrix
         cout << "Building new Density matrix, Pnew" << endl;
         //Matrix P = Matrix::Zero(ao,ao);
         Matrix C_ao_new = Matrix::Zero(ao,ao);
	 Matrix P = Matrix::Zero(ao,ao);
         diagonalize_Fock(ao, Fock_new, Xmat, C_ao_new, evals, evecs);
         build_P(ao, occ, C_ao_new, P);
	 //cout << "P after function: \n" << P << endl;
      //Compute new SCF Energy
         En_elec_new =calculate_En_elec(ao, P, H_core, Fock_new);
         En_total = En_elec_new + En_nuc;
	 //cout << "P after energy: \n" << P << endl;
         cout << "Total Energy is ...." << endl;
         cout << En_elec_new << "+" << En_nuc << "=" << En_total << endl;
	 
	 iteration = iteration + 1;
      
         //cout << "The energy is: " << hf_energy << endl;
         deltaE = En_elec_new-En_elec;
	 P2=P;
         cout << "deltaE is: " << deltaE << endl;
	 C_mo = C_ao_new;
	 cout << "C_ao_new at end of while\n" << C_ao_new << endl;
   }

//Begin MP2 here
   cout << "C_ao_new at beginning\n" << C_mo << endl;
   cout << "Beginning MP2" <<  endl;

//transformation
   cout << "Transforming to MO basis" << endl;
   Real_4dMatrix v_int_mo(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
   Real_4dMatrix v_int_mo_2(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
   transform_v_int(ao, C_mo, v_int, v_int_mo, Xmat, evecs);
   transform_v_int_2(ao, v_int, v_int_mo_2, Xmat, C_mo);
//calculate MP2 energy
   cout << "------------------------------------Calculating MP2 energy------------------------------------" << endl;
////////we want to get Emp2 = -0.049149636120
   double Emp2 = calculate_E_mp2(ao, occ, evals, v_int_mo);
   double Emp2_2 = calculate_E_mp2(ao, occ, evals, v_int_mo_2);

   cout << "The final energy is: " << En_elec_new + Emp2 + En_nuc <<endl; 
   cout << "The final energy with the N^5 v_int_mo is: " << En_elec_new + Emp2_2 + En_nuc << endl;
   return 0;
}


