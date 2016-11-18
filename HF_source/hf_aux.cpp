#include <iostream>
#include <fstream>
#include <iomanip>
//#include <cstdio>
#include "hf_aux.hpp"
#include <string>
#include <string.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
//////Here are the functions required for HF
using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;

double read_nuc_en(char* Path)
{
   string ss(Path);
   string name("");
   name=ss +"/enuc.dat";
   double nuc_en;
   std::ifstream nuc_ener;
   nuc_ener.open(name.c_str());
   if (!nuc_ener.is_open()) {
     cout << "Error: input file nuc_ener cannot open"<< endl;
   }
   while(!nuc_ener.eof()){
      nuc_ener >> nuc_en;
   }
  // cout << "in read_nuc routine en = " << nuc_en << endl;

   return nuc_en;

}


void read_T(int ao, Matrix& T_int,char* Path){
   
   double val;
   int i;
   int j;
   
   string ss(Path);
   string name("");
   name=ss +"/T.dat";

   std::ifstream kin_en;
   kin_en.open(name.c_str());
   if (!kin_en.is_open()) {
     cout << "Error: input file kin_en cannot open"<< endl;
   }
   for(int k=0;k<((ao)*(ao+1))/2;k++){
      kin_en >> i;
      kin_en >> j;
      kin_en >> val;
      T_int(i-1,j-1) = val;
      T_int(j-1,i-1) = T_int(i-1,j-1);
      //std::cout << i << " " << j << " " << T_int(i-1,j-1) << std::endl;
   }
   return;
}

void read_S(int ao, Matrix& S,char* Path){

   double val;
   int i;
   int j;

   string ss(Path);              
   string name("");    
   name=ss +"/overlap.dat";

   std::ifstream overlap;
   overlap.open(name.c_str());
   if (!overlap.is_open()) {
     cout << "Error: input file overlap cannot open"<< endl;
   }
   for(int k=0;k<((ao)*(ao+1))/2;k++){
      overlap >> i;
      overlap >> j;
      overlap >> val;
      S(i-1,j-1) = val;
      S(j-1,i-1) = S(i-1,j-1);
      //std::cout << i << " " << j << " " << S(i-1,j-1) << std::endl;
   }
   return;
}

void read_v_nuc(int ao, Matrix& v_nuc, char* Path){

   double val;
   int i;
   int j;

   string ss(Path);     
   string name("");  
   name=ss +"/v_nuc.dat";

   std::ifstream interaction;
   interaction.open(name.c_str());
   if (!interaction.is_open()) {
     cout << "Error: input file v_nuc cannot open"<< endl;
   }
   for(int k=0;k<((ao)*(ao+1))/2;k++){
      interaction >> i;
      interaction >> j;
      interaction >> val;
      v_nuc(i-1,j-1) = val;
      v_nuc(j-1,i-1) = v_nuc(i-1,j-1);
      //std::cout << i << " " << j << " " << v_nuc(i-1,j-1) << std::endl;
   }
   return;

}


void build_H_core(int ao, Matrix& v_nuc, Matrix& T_int, Matrix& H_core){

   H_core=T_int+v_nuc;
   //std::cout << H_core << endl;

   return;
}

void read_v_int(int ao, Matrix& V, char* Path){

   double val;
   int i;
   int j;
   int k;
   int l;
   int z;
  
   string ss(Path);                      
   string name("");        
   name=ss +"/v_int.dat";

   std::ifstream interaction;
   interaction.open(name.c_str());
   if (!interaction.is_open()) {
     cout << "Error: input file interaction cannot open"<< endl;
   }
   while(!interaction.eof()){      
      interaction >> i;
      interaction >> j;
      interaction >> k;
      interaction >> l;
      interaction >> val;
      V((i-1)*ao+j-1,(k-1)*ao+l-1)=val;



   }

   for (int mu=0; mu<ao;mu++){
      for (int nu=0; nu<mu+1;nu++){
          for (int lam=0; lam<ao;lam++){
             for(int sig=0; sig<lam+1;sig++){
                   if (mu*(mu+1)/2 + nu >=lam*(lam+1)/2+sig){
                        V(nu*ao+mu,lam*ao+sig)=V(mu*ao+nu,lam*ao+sig);
                        V(mu*ao+nu,sig*ao+lam)=V(mu*ao+nu,lam*ao+sig);
                        V(lam*ao+sig,mu*ao+nu)=V(mu*ao+nu,lam*ao+sig);
                        V(sig*ao+lam,mu*ao+nu)=V(mu*ao+nu,lam*ao+sig);
                        V(lam*ao+sig,nu*ao+mu)=V(mu*ao+nu,lam*ao+sig);
                        V(sig*ao+lam,nu*ao+mu)=V(mu*ao+nu,lam*ao+sig);
                        V(nu*ao+mu,sig*ao+lam)=V(mu*ao+nu,lam*ao+sig); 

                    }


             }
         }
      }
   }

   return;
};

void calculate_S12(int ao, Matrix& S, Matrix& S12, Matrix& Xmat){

   Eigen::SelfAdjointEigenSolver<Matrix> solver(S);
   Matrix evecs = solver.eigenvectors();
   Matrix evals = solver.eigenvalues();

   //std::cout << "eigenvectors are: " << evecs << std::endl;
   //std::cout << "eigenvalues are: " << evals << std::endl;

   for (int i=0; i < ao ; i++){
      //for (int j=0; j < ao ; j++){
          //std::cout << evals(i) << std::endl;
          S12(i,i)=1/sqrt(evals(i));
      //}
   }
   
   //std::cout << "Square root of eigenvalues are: " << S12 << std::endl;
   
   //Matrix evecs_trans = evecs.transpose();
   Xmat = evecs*S12*evecs.transpose();

   //std::cout << "Xmat is: " << Xmat << endl;


   return;
};


void diagonalize_Fock(int ao, Matrix& H_core, Matrix& Xmat, Matrix& C_ao, Matrix& evals, Matrix& evecs){

   Matrix Fock = Matrix::Zero(ao,ao);
   Fock=Xmat.transpose()*H_core*Xmat;
      
   Eigen::SelfAdjointEigenSolver<Matrix> solver(Fock);
   evecs = solver.eigenvectors();
   evals = solver.eigenvalues();

   //cout << "this is what evals looks like" << evals << endl;
   //cout << "this is what evecs looks like" << evecs << endl;
   C_ao=Xmat*evecs;

   //cout << "C_ao in AO basis is: " << C_ao << endl;


   return;

};



void build_P(int ao, int occ, Matrix& C_ao, Matrix &P0){

  for (int i=0; i<ao; i++) {
      for (int j=0; j<ao; j++) {
         for (int k=0; k<occ; k++) {
       
            P0(i,j) = P0(i,j) + C_ao(i,k)*(C_ao(j,k));
         
         }
      }
   }

   //cout << "Initial Density matrix is: " << P0 << endl; 
   return;
};

double calculate_En_elec(int ao, Matrix& P0, Matrix& H_core, Matrix& Fock){
  double En = 0.0;
  for (int i=0; i<ao; i++) {
      for (int j=0; j<ao; j++) {
      
          En = En + P0(j,i)*(H_core(i,j)+Fock(i,j));

      }
   }

   //cout << "En_elec in loop is: " << En << endl;
   return En;
};


void build_new_Fock(int ao, Matrix& P0, Matrix& v_int, Matrix&H_core, Matrix&Fock){

   for (int i=0; i<ao ; i++){
      for (int j=0; j<ao ; j++){
         Fock(i,j)=H_core(i,j);
         for (int k=0; k<ao ; k++){
             for (int l=0; l<ao ; l++){

                 Fock(i,j)=Fock(i,j) +  P0(k,l)*(2*v_int(i*ao+j,k*ao+l)-v_int(i*ao+k,j*ao+l));

             }
          }
      }
   }

   //cout << "New fock (G) = " << endl << Fock << endl;

};

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
	 C_mo = C_ao_new;

   }
   return En_total;
}

