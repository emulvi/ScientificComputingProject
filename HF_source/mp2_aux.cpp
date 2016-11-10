#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>



typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;

void transform_v_int(int ao, Matrix& C, Real_4dMatrix& v_int, Real_4dMatrix& v_int_mo, Matrix& Xmat, Matrix& evecs){

   int i, j, k, l;
   int a,b,c,d = 0;
 
   //C=C.transpose();
 
   //Matrix Cmo = Xmat.transpose()*C;
   Matrix Cmo = evecs.transpose();
   //cout << "C matrix is: " << Cmo << endl;
   cout << "helllooo? " << endl;
 
   for(a=0; a < ao; a++) {
     for(b=0; b < ao; b++) {
       for(c=0; c < ao; c++) {
         for(d=0; d < ao; d++) {
  
           for(i=0; i < ao; i++) {
             for(j=0; j < ao; j++) {
               for(k=0; k < ao; k++) {
                 for(l=0; l < ao; l++) {
  
                   v_int_mo[a][b][c][d] += Cmo(i,a) * Cmo(j,b) * Cmo(k,c) * Cmo(l,d) * v_int[i][j][k][l];
                 }
               }
             }
           }
  
         }
       }
     }
   }

  
   return;
};


double calculate_E_mp2(int ao, int occ, Matrix& evals, Real_4dMatrix& v_int_mo){
   double E_mp2 = 0.0;

   for(int i=0; i < occ; i++) {
     for(int a=occ; a < ao; a++) {
       for(int j=0; j < occ; j++) {
         for(int b=occ; b < ao; b++) {

            E_mp2 += v_int_mo[i][a][j][b]*((2*v_int_mo[i][a][j][b]-v_int_mo[i][b][j][a])/(evals(i)+evals(j)-evals(a)-evals(b)));

         }
       }
     }
   }

   cout << "MP2 energy: " << E_mp2 << endl;
   return E_mp2;
};
