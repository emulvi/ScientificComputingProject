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
//int nthreads;
//int tid;
int threads = 8;
double Pi_seed[threads];

//#pragma omp parallel private(nthreads,tid)
#pragma omp parallel for
for (int i=0; i<threads; ++i){
    //tid = omp_get_thread_num();
    Pi_seed[i]=findPi(N);
    cout << "Pi at processor " << i << " is: " <<Pi_seed[i] << endl;
  /* Only master thread does this */
   // if (tid == 0) 
    //{
    //nthreads = omp_get_num_threads();
    //printf("Number of threads = %d\n", nthreads);

    //}

    }
for (int i=0; i<threads; ++i){
cout << "Pi at processor " << i << " is: " << Pi_seed[i] << endl;
}
double sum_pi=0;
for (int i=0; i<threads; ++i){
    sum_pi = sum_pi + Pi_seed[i];
}

cout << "The final value of pi is " << sum_pi /threads << endl;


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




