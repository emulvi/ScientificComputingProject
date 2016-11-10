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

void transform_v_int(int ao, Matrix& C, Real_4dMatrix& v_int, Real_4dMatrix& v_int_mo){

  int i, j, k, l;
  int a,b,c,d;

  for(a=0; a < ao; a++) {
    for(b=0; b <= a; b++) {
      for(c=0; c <= a; c++) {
        for(d=0; d <= (a==c ? b : c); d++) {
 
          for(i=0; i < ao; i++) {
            for(j=0; j < ao; j++) {
              for(k=0; k < ao; k++) {
                for(l=0; l < ao; l++) {
 
                  v_int_mo[i][j][k][l] += C(d,l) * C(c,i) * C(b,j) * C(a,i) * v_int[a][b][c][d];
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

   for(int a=0; a < occ; a++) {
     for(int b=0; b < occ; b++) {
       for(int s=occ; s < ao; s++) {
         for(int r=occ; r < ao; r++) {

    E_mp2 += v_int_mo[a][s][b][r]*((2*v_int_mo[a][s][b][r]-v_int_mo[a][r][b][s])/(evals(a)+evals(b)-evals(r)-evals(s)));

         }
       }
     }
   }

   cout << "MP2 energy: " << E_mp2 << endl;
   return E_mp2;
};
