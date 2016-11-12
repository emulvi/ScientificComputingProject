#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "Eigen/Cholesky"
#include <iomanip>



typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;

void transform_v_int(int ao, Matrix& C, Matrix& v_int, Matrix& v_int_mo, Matrix& Xmat, Matrix& evecs){

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
  
                   v_int_mo(a*ao+b,c*ao+d) += C(i,a) * C(j,b) * C(k,c) * C(l,d)*v_int(i*ao+j,k*ao+l); 
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

void transform_v_int_2(int ao, Matrix& v_int, Matrix& v_int_mo_2, Matrix& Xmat, Matrix& C_mo){
  int a,b,c,d,e;
  //Real_4dMatrix v1(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
  //Real_4dMatrix v2(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
  //Real_4dMatrix v3(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
  Matrix v1 = Matrix::Zero(ao*ao,ao*ao);
  Matrix v2 = Matrix::Zero(ao*ao,ao*ao);
  Matrix v3 = Matrix::Zero(ao*ao,ao*ao);

  std::clock_t start;
  start = std::clock();

  for (a=0;a<ao;a++){
    for (b=0;b<ao;b++){
      for (c=0;c<ao;c++){	
	for (d=0;d<ao;d++){
	  for (e=0;e<ao;e++){
	    v1(a*ao+b,c*ao+d) += v_int(a*ao+b,c*ao+e)*C_mo(e,d);
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
	    v2(a*ao+b,c*ao+d) += v1(a*ao+b,e*ao+d)*C_mo(e,c);
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
	    v3(a*ao+b,c*ao+d) += v2(a*ao+e,c*ao+d)*C_mo(e,b);
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
	    v_int_mo_2(a*ao+b,c*ao+d) += v3(e*ao+b,c*ao+d)*C_mo(e,a);
	  }
	}
      }
    }
  }

   std::cout << "Time in N^5: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;

  return;
}


void transform_v_int_CD(int ao, Matrix& v_int, Matrix& v_int_mo_2, Matrix& Xmat, Matrix& C_mo){
  int a,b,c,d,e;
//  Real_4dMatrix v1(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
//  Real_4dMatrix v2(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
//  Real_4dMatrix v3(ao, vector<vector<vector<double> > >(ao, vector<vector<double> >(ao, vector<double>(ao,0.0))));
//
//  std::clock_t start;
//  start = std::clock();
//
//
//  std::cout << "Cholesky Decomp" << endl;
//
//  //Eigen::LDLT< Real_Matrix, Upper > ldlt_object(C_mo);
//  //Eigen::Matrix C_mo(ao,ao);
//  Eigen::LDLT<Matrix> ldltOfC_mo(C_mo);
//  Matrix L = ldltOfC_mo.matrixL();
//  Matrix D = ldltOfC_mo.vectorD();
//  //Eigen::Diagonal<double Eigen::MatrixXd>vectorD() double;
//  //Eigen::VectorXd D(C_mo.llt().vectorD() );
//
//  std::cout << "L is: " << L << std::endl;
//  std::cout << "L^T is: " << L.transpose() << std::endl;
//  std::cout << "D is: " << D << std::endl;
//
//  std::cout << "Time in N^4: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;

  return;
}


double calculate_E_mp2(int ao, int occ, Matrix& evals, Matrix& v_int_mo){
   double E_mp2 = 0.0;

   //cout << "evals in mp2 are: " << evals << endl;
   int docc=5;
   occ=docc;

   for(int i=0; i < occ; i++) {
       for(int j=0; j < occ; j++) {
	 for(int a=occ; a < ao; a++) {
	   for(int b=occ; b < ao; b++) {

            E_mp2 += v_int_mo(i*ao+a,j*ao+b)*(2*v_int_mo(i*ao+a,j*ao+b)-v_int_mo(i*ao+b,j*ao+a))/(evals(i)+evals(j)-evals(a)-evals(b));

         }
       }
     }
   }

   cout << "MP2 energy: " << E_mp2 << endl;
   return E_mp2;
};
