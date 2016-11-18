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
//#include "hf_main.cpp"

////note, must include -Wno-write-strings in test command
using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;


void test_En_nuc();
void test_HF_en();




int main()
{

test_En_nuc();
test_HF_en();


}



void test_HF_en(){

int ao = 7;
int occ = 5;
char *path ="Water_STO-3G";

double En = hf_main(ao,occ,path);

if (abs(abs(En)-74.9421)<=pow(10,-4)){cout << "HF_Energy:  ... Test passed" << endl;}else{cout << "HF_Energy:  ... Test failed" << endl;}

}



void test_En_nuc(){
int ao = 7;
int occ = 5;
char *path = "Water_STO-3G";


double En_nuc = read_nuc_en(path);

if (abs(abs(En_nuc)-8.002367061810450)<=pow(10,-8)){cout << "Nuclear_Energy:  ...Test passed" << endl;}else{cout << "Nuclear_Energy:  ...Test failed" << endl;}

};
