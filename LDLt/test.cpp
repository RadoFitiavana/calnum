#include "tools.h"
#include <cstdlib>
#include<iostream>
using namespace std ;

int main(){
	char readkey= 'n' ;
	while (readkey == 'n'){
	    system("cls") ;
	    Factor f ;
	    f.factorize() ;
	    f.solve() ;
	    f.display() ;
	    cout <<"Close [y|n]: "; cin>>readkey;
    }
    return 0 ;
}