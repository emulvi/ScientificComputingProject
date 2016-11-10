#include <iostream>
#include <fstream>
#include <iomanip>
//#include <cstdio>
#include "hf_aux.hpp"
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
//////Here are the functions required for HF


typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;

double read_nuc_en()
{
   double nuc_en;
   std::ifstream nuc_ener("enuc.dat");
   while(!nuc_ener.eof()){
      nuc_ener >> nuc_en;
   }
   cout << "in read_nuc routine en = " << nuc_en << endl;

   return nuc_en;

}


void read_T(int ao, Matrix& T_int){
   
   double val;
   int i;
   int j;

   std::ifstream kin_en;
   kin_en.open("T.dat");
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

void read_S(int ao, Matrix& S){

   double val;
   int i;
   int j;

   std::ifstream overlap;
   overlap.open("overlap.dat");
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

void read_v_nuc(int ao, Matrix& v_nuc){

   double val;
   int i;
   int j;


   std::ifstream interaction;
   interaction.open("v_nuc.dat");
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
   std::cout << H_core << endl;

   return;
}

void read_v_int(int ao, Real_4dMatrix& v_int){

   double val;
   int i;
   int j;
   int k;
   int l;
   int z;
  

   std::ifstream interaction;
   interaction.open("v_int.dat");
   while(!interaction.eof()){      
      interaction >> i;
      interaction >> j;
      interaction >> k;
      interaction >> l;
      interaction >> val;
      v_int[i-1][j-1][k-1][l-1] = val;
      //std::cout << i << " " << j << " " << k << " " << l << " " << v_int[i-1][j-1][k-1][l-1] << std::endl;
   }

   for (int mu=0; mu<ao;mu++){
      for (int nu=0; nu<mu+1;nu++){
          for (int lam=0; lam<ao;lam++){
             for(int sig=0; sig<lam+1;sig++){
                   if (mu*(mu+1)/2 + nu >=lam*(lam+1)/2+sig){
                        v_int[nu][mu][lam][sig]=v_int[mu][nu][lam][sig];
                        v_int[mu][nu][sig][lam]=v_int[mu][nu][lam][sig];
                        v_int[nu][mu][sig][lam]=v_int[mu][nu][lam][sig];
                        //v_int[nu][mu][lam][sig]=v_int[mu][nu][lam][sig];
                        v_int[lam][sig][mu][nu]=v_int[mu][nu][lam][sig];
                        v_int[sig][lam][mu][nu]=v_int[mu][nu][lam][sig];
                        v_int[lam][sig][nu][mu]=v_int[mu][nu][lam][sig];
                        v_int[sig][lam][nu][mu]=v_int[mu][nu][lam][sig];
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
          std::cout << evals(i) << std::endl;
          S12(i,i)=1/sqrt(evals(i));
      //}
   }
   
   //std::cout << "Square root of eigenvalues are: " << S12 << std::endl;
   
   //Matrix evecs_trans = evecs.transpose();
   Xmat = evecs*S12*evecs.transpose();

   //std::cout << "Xmat is: " << Xmat << endl;


   return;
};


void diagonalize_Fock(int ao, Matrix& H_core, Matrix& Xmat, Matrix& C_ao){

   Matrix Fock = Matrix::Zero(ao,ao);
   Fock=Xmat.transpose()*H_core*Xmat;
      
   Eigen::SelfAdjointEigenSolver<Matrix> solver(Fock);
   Matrix evecs = solver.eigenvectors();
   Matrix evals = solver.eigenvalues();

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

   cout << "Initial Density matrix is: " << P0 << endl; 
   return;
};

double calculate_En_elec(int ao, Matrix& P0, Matrix& H_core, Matrix& Fock){
  double En = 0.0;
  for (int i=0; i<ao; i++) {
      for (int j=0; j<ao; j++) {
      
          En = En + P0(j,i)*(H_core(i,j)+Fock(i,j));

      }
   }

   cout << "En_elec in loop is: " << En << endl;
   return En;
};


void build_new_Fock(int ao, Matrix& P0, Real_4dMatrix& v_int, Matrix&H_core, Matrix&Fock){

   for (int i=0; i<ao ; i++){
      for (int j=0; j<ao ; j++){
         Fock(i,j)=H_core(i,j);
         for (int k=0; k<ao ; k++){
             for (int l=0; l<ao ; l++){

                 Fock(i,j)=Fock(i,j) +  P0(k,l)*(2*v_int[i][j][k][l]-v_int[i][k][j][l]);

             }
          }
      }
   }

   //cout << "New fock (G) = " << endl << Fock << endl;

}

