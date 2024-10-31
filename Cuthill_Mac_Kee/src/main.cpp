#include "tools.h"
#include <iostream>
#include <cstdlib>

using namespace std ;

void call_Solver() ;

int main(){
	system("clear") ;
	call_Solver() ;
	return 0 ;
}

void call_Solver(){
	Solver S ;
	S.factorize() ;
	S.solve() ;
	S.display() ;
}