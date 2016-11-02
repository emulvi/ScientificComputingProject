#include <iostream>
#include <time.h> // time
#include <stdlib.h> // srand
#include <iomanip>
#include <omp.h> 
using namespace std;

int main(int argc, char *argv[]){

double findPi(const int);
double Pi;
const int N={1e5};
int nthreads;
int tid;
double Pi_seed[nthreads];

#pragma omp parallel for(nthreads,tid)

    {
    tid = omp_get_thread_num();
    Pi_seed[tid]=findPi(N);
    cout << "Pi at processor " << tid << " is: " <<Pi_seed[tid] << endl;}

}
double findPi(const int N) {
    srand(time(NULL));
    cout.precision(100);
    srand(time(0));
    int circle = 0;
    for (int i = 0; i < N; i ++) {
        double x = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        double y = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        if (x * x + y * y <= 1.0) circle ++;
    }

    double Pi=(double)circle / N * 4.0;
//    cout << N << (char)9 << (char)9 << Pi << endl;
    return Pi;
}




