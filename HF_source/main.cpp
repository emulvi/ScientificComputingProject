#include <iostream>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <string.h>
#include <stdio.h>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Cholesky"
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

   for (int i=1; i<7; i=i+2){

     if (strncmp(argv[i],"-ao",3)==0) {
         ao =  atoi(argv[i + 1]);
     } 
     else if (strncmp(argv[i],"-occ",4)==0) {
       occ = atoi(argv[i + 1]);
     } 
     else if(strncmp(argv[i], "-path",5)==0) {
       Path = argv[i + 1];
     }
     else {
       cout << "The most common command-line options are: '-ao 7 -occ 5 -path Water_STO-3G'\n Try again with new input parameters......\n Now exiting. " << endl;
       exit(0);
     };
   };
   
   cout << "Initial Values" << endl;
   cout << "\t" << "ao:     " << ao << endl;
   cout << "\t" << "occ:    " << occ << endl;
   cout << "\t" << "Path:   " << Path << endl;
   
   Matrix C_mo = Matrix::Zero(ao,ao);
   Matrix V = Matrix::Zero(ao*ao,ao*ao);
   Matrix evals = Matrix::Zero(ao,ao);
   double En_total = 0;
   double Emp2 = 0;

   En_total = hf_main(ao, occ, Path, V, C_mo, evals);

   cout << "Beginning MP2" <<  endl;
   
   Emp2 = mp2_main(ao, occ, Path, V, C_mo, evals, En_total);

   return 0;
}


