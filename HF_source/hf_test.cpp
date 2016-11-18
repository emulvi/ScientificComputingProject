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

////note, must include -Wno-write-strings in test command
using namespace std;
//using namespace Eigen;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef std::vector<vector<double> > Real_Matrix;
typedef std::vector<vector<vector<vector<double> > > > Real_4dMatrix;


int main()
{


int ao = 7;
int occ = 5;
char *path = "Water_STO-3G";


double En_nuc = read_nuc_en(path);
std::cout << "Nuclear energy is: " << En_nuc << std::endl;

if (En_nuc-8.002367061810450<=pow(10,-8)){cout << "Test passed" << endl;}


}
