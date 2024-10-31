#include <iostream>
#include "range_profil.h"
#include <cstdlib>
using namespace std ;

int main(){
	char readkey('n') ;
	while(readkey == 'n'){
	    system("cls") ;
	    Profil P ;
	    P.factorize() ;
	    P.solve() ;
	    P.display() ;
	    cout<<"Close [y|n]: " ; cin>>readkey ;
	} 
	return 0 ;
}