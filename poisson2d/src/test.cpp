#include "poisson.h"
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std ;

int main() {
    cout << Poisson::nCores << endl ;
    Poisson::ParallSetup::setPARALLEL(true) ;
    srand(time(NULL)) ;
    int N = 15 ;
    try{
        float *sgn = new float[N] ;
        for (int i=0; i<N; i++) {
            sgn[i] = (rand()%2)+1 ;
            cout << sgn[i] <<" " ;
        }
        Poisson::FDST<float> fst (sgn,N) ;
        fst.runFDST() ;
        cout << "\n" <<endl ;
        for (int i=0; i<N; i++) {
            cout << sgn[i] << " " ;
        }
        cout <<endl ;
    }
    catch(exception& e) {
        cout << e.what() << endl ;
    }
    return 0 ;
}
