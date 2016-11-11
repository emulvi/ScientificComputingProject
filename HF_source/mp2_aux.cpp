#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>



typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;

void transform_v_int(int ao, Matrix& C, Real_4dMatrix& v_int, Real_4dMatrix& v_int_mo, Matrix& Xmat, Matrix& evecs){

   int i, j, k, l;
   int a,b,c,d;
 
   std::clock_t start;
   start = std::clock();
 
   for(a=0; a < ao; a++) {
     for(b=0; b < ao; b++) {
       for(c=0; c < ao; c++) {
         for(d=0; d < ao; d++) {
  
           for(i=0; i < ao; i++) {
             for(j=0; j < ao; j++) {
               for(k=0; k < ao; k++) {
                 for(l=0; l < ao; l++) {
  
                   v_int_mo[a][b][c][d] += C(i,a) * C(j,b) * C(k,c) * C(l,d)*v_int[i][j][k][l]; 
		   //cout << v_int_mo[a][b][b][d] << endl;
                 }
               }
             }
           }  
         }
       }
     }
   }

   std::cout << "Time in N^8: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;  
   return;
};

void transform_v_int_2(int ao, Real_4dMatrix& v_int, Real_4dMatrix& v_int_mo_2, Matrix& Xmat, Matrix& C_mo){
  int a,b,c,d,e;
  Real_4dMatrix v1(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
  Real_4dMatrix v2(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
  Real_4dMatrix v3(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));

  std::clock_t start;
  start = std::clock();

  for (a=0;a<ao;a++){
    for (b=0;b<ao;b++){
      for (c=0;c<ao;c++){	
	for (d=0;d<ao;d++){
	  for (e=0;e<ao;e++){
	    v1[a][b][c][d] += v_int[a][b][c][e]*C_mo(e,d);
	  }
	}
      }
    }
  }
  for (a=0;a<ao;a++){
    for (b=0;b<ao;b++){
      for (c=0;c<ao;c++){
	for (d=0;d<ao;d++){
	  for (e=0;e<ao;e++){
	    v2[a][b][c][d] += v1[a][b][e][d]*C_mo(e,c);
	  }
	}
      }
    }
  }
  for (a=0;a<ao;a++){
    for (b=0;b<ao;b++){
      for (c=0;c<ao;c++){
	for (d=0;d<ao;d++){
	  for (e=0;e<ao;e++){
	    v3[a][b][c][d] += v2[a][e][c][d]*C_mo(e,b);
	  }
	}
      }
    }
  }
  for (a=0;a<ao;a++){
    for (b=0;b<ao;b++){
      for (c=0;c<ao;c++){
	for (d=0;d<ao;d++){
	  for (e=0;e<ao;e++){
	    v_int_mo_2[a][b][c][d] += v3[e][b][c][d]*C_mo(e,a);
	  }
	}
      }
    }
  }

   std::cout << "Time in N^5: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;

  return;
}


double calculate_E_mp2(int ao, int occ, Matrix& evals, Real_4dMatrix& v_int_mo){
   double E_mp2 = 0.0;

   //cout << "evals in mp2 are: " << evals << endl;
   int docc=5;
   occ=docc;

   for(int i=0; i < occ; i++) {
       for(int j=0; j < occ; j++) {
	 for(int a=occ; a < ao; a++) {
	   for(int b=occ; b < ao; b++) {

            E_mp2 += v_int_mo[i][a][j][b]*(2*v_int_mo[i][a][j][b]-v_int_mo[i][b][j][a])/(evals(i)+evals(j)-evals(a)-evals(b));

         }
       }
     }
   }

   cout << "MP2 energy: " << E_mp2 << endl;
   return E_mp2;
};
