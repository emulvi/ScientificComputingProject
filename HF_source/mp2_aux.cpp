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

void transform_v_int_3(int occ, int ao, Matrix& v_int, Matrix& v_int_mo_2, Matrix& Xmat, Matrix& C_mo){
  int a,b,c,d,e;


  Matrix v1 = Matrix::Zero(ao*ao,ao*ao);
  Matrix v2 = Matrix::Zero(ao*ao,ao*ao);
  Matrix v3 = Matrix::Zero(ao*ao,ao*ao);

  std::clock_t start;
  start = std::clock();


  for (a=0;a<ao;a++){
    for (b=0;b<ao;b++){
      for (c=0;c<ao;c++){	
	for (d=occ;d<ao;d++){
	  for (e=0;e<ao;e++){
	    v1(a*ao+b,c*ao+d) += v_int(a*ao+b,c*ao+e)*C_mo(e,d);
	  }
	}
      }
    }
  }
  for (a=0;a<ao;a++){
    for (b=0;b<ao;b++){
      for (c=0;c<occ;c++){
	for (d=occ;d<ao;d++){
	  for (e=0;e<ao;e++){
	    v2(a*ao+b,c*ao+d) += v1(a*ao+b,e*ao+d)*C_mo(e,c);
	  }
	}
      }
    }
  }
  for (a=0;a<ao;a++){
    for (b=occ;b<ao;b++){
      for (c=0;c<occ;c++){
	for (d=occ;d<ao;d++){
	  for (e=0;e<ao;e++){
	    v3(a*ao+b,c*ao+d) += v2(a*ao+e,c*ao+d)*C_mo(e,b);
	  }
	}
      }
    }
  }
  for (a=0;a<occ;a++){
    for (b=occ;b<ao;b++){
      for (c=0;c<occ;c++){
	for (d=occ;d<ao;d++){
	  for (e=0;e<ao;e++){
	    v_int_mo_2(a*ao+b,c*ao+d) += v3(e*ao+b,c*ao+d)*C_mo(e,a);
	  }
	}
      }
    }
  }

   std::cout << "Time in v^2o^2: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;

  return;
}


void transform_v_int_CD(int ao, Matrix& v_int, Matrix& v_int_mo_2, Matrix& Xmat, Matrix& C_mo){
  int a,b,c,d,e;

  Matrix v1 = Matrix::Zero(ao*ao,ao*ao);
  Matrix v2 = Matrix::Zero(ao*ao,ao*ao);
  Matrix v3 = Matrix::Zero(ao*ao,ao*ao);

  std::clock_t start;
  start = std::clock();


  std::cout << "Cholesky Decomp" << endl;

  Eigen::LLT<Matrix> lltOfC_mo(C_mo);
  Matrix L = lltOfC_mo.matrixL();
//  Matrix D = ldltOfC_mo.vectorD();

  std::cout << "L is: " << L << std::endl;
  std::cout << "L^T is: " << L.transpose() << std::endl;
  std::cout << "the shape of L is: " << L.size() << L.cols() << L.rows() << std::endl;

//  std::cout << "D is: " << D << std::endl;

//  Eigen::LLT<Matrix> lltOfL(L);
//  Matrix LL = lltOfL.matrixL();
//
//  Eigen::LLT<Matrix> lltOfLT(L.transpose());
//  Matrix LT = lltOfLT.matrixL();
////  Matrix D = ldltOfC_mo.vectorD();
//
//
//  std::cout << "LL is: " << LL << std::endl;
//  std::cout << "LL^T is: " << LL.transpose() << std::endl;
//
//  std::cout << "LT is: " << LT << std::endl;
//  std::cout << "LT^T is: " << LT.transpose() << std::endl;

  Matrix tempA = Matrix::Zero(ao,ao);
  Matrix tempB = Matrix::Zero(ao,ao);
  Matrix tmpLT = L.transpose();
  for (int i=0;i<ao;i++){
    for (int j=0;j<ao;j++){
      for (int sig=0;sig<ao;sig++){
        for (int rho=0;rho<ao;rho++){
           tempA(i,j) += C_mo(sig,i)*C_mo(rho,j)*tmpLT(sig,rho);
        }
      }
    }
  }
  

  for (int a=0;a<ao;a++){
    for (int b=0;b<ao;b++){
      for (int mu=0;mu<ao;mu++){
        for (int nu=0;nu<ao;nu++){
           tempB(a,b) += C_mo(a,mu)*C_mo(b,nu)*L(mu,nu)*tempA(a,b);
        }
      }
    }
  }

  for (int a=0;a<ao;a++){
    for (int b=0;b<ao;b++){
      for (int i=0;i<ao;i++){
        for (int j=0;j<ao;j++){
           v_int_mo_2(a*ao+b,i*ao+j) += tempB(a,b)*tempA(i,j);
        }
      }
    }
  }

  std::cout << "Time in N^4: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;
  return;
}


double calculate_E_mp2(int ao, int occ, Matrix& evals, Matrix& v_int_mo){
   double E_mp2 = 0.0;

   //cout << "evals in mp2 are: " << evals << endl;
   int docc=5;
   occ=docc;

   cout << "hello" << endl;
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
