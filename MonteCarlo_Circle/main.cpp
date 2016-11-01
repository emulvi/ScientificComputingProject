#include <iostream>
#include <time.h> // time
#include <stdlib.h> // srand
 
using namespace std;

double findPi(const int);

int main(){
const int N={1e3};
int num_proc=5;
double Pi_seed[num_proc];
    for (int j=0; j<num_proc; j++){
        Pi_seed[j]=findPi(N);
        cout << "Pi at processor " << j << " is: " <<Pi_seed[j] << endl;
    }
}


double findPi(const int N) {
    srand(time(NULL));
    cout.precision(10);
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




