#include <iostream>

using namespace std ;

int main() {
	char *ptr = new char ;
	char *p = ptr ;
	*ptr = 'a' ;
	cout << "initial ptr: "<< *p << endl ;
	delete ptr ;
	*(ptr+1) = 'b' ;
	p = ptr+1 ;
	cout << "current ptr+1: " << *p << endl ;
	//*p = 'c' ;
	//cout << "modified ptr: " << *ptr << endl ;
	return 0 ;
}
