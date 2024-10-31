#include <iostream>
#include <cmath>
#include <fstream>
#include "tools.h"
using namespace std ;

Factor::Factor(){
	ifstream f("matrix.txt") ;
	if (f){
		f >> dim ;
		try{
			A = new double*[dim] ;
			L = new double*[dim] ;
			D = new double[dim] ;
			b = new double[dim] ;
			x = new double[dim] ;
			for (int i(0); i<dim; i++){
				A[i] = new double[i+1] ;
				L[i] = new double[i+1] ;
				for (int j(0); j<i+1; j++){
					f >> A[i][j] ;
				}
			}
			for (int i(0); i<dim; i++){
				f >> b[i] ;
			}
		}
		catch(exception& e){
			cout << e.what()<<endl ;
		}
		f.close() ;
	}
	else{
		cout << "Error while opening file matrix.txt"<<endl;
	}
}

Factor::~Factor(){
	for (int i(0); i<dim; i++){
		delete[] A[i] ;
		delete[] L[i] ;
	}
	delete[] b ;
	delete[] x ;
	delete[] D ;
	delete[] A ;
	delete[] L ; 
}

void Factor::factorize(){
	D[0] = A[0][0] ;
	double r ;
	L[0][0] = 1.0 ;
	for (int i(1); i<dim; i++){
		L[i][i] = 1 ;
		for (int j(0); j<i; j++){
			 r = 0.0 ;
			for (int k(0); k<j; k++){
				r += L[i][k]*D[k]*L[j][k] ;
			}
			L[i][j] = A[i][j] - r ;
			L[i][j] /= D[j] ;
		}
		r = 0.0 ;
		for (int k(0); k<i; k++){
			r += D[k]*L[i][k]*L[i][k] ;
		}
		D[i] = A[i][i] - r ;
	}
}

void Factor::display(){
	cout << "Matrix A:\n-----------"<<endl ;
	for (int i(0); i<dim; i++){
		for (int j(0); j<i+1; j++){
			cout << A[i][j]<<" " ;
		}
		for (int j(i+1); j<dim; j++){
			cout << A[j][i]<<" " ;
		}
		cout <<endl ;
	}
	cout <<endl ;	
	cout << "Matrix L:\n-----------"<<endl ;
	for (int i(0); i<dim; i++){
		for (int j(0); j<i+1; j++){
			cout << L[i][j]<<" " ;
		}
		for (int j(i+1); j<dim; j++){
			cout << L[j][i]<<" " ;
		}
		cout <<endl ;
	}
	cout <<endl ;
	cout << "Matrix D:\n-----------"<<endl ;
	for (int i(0); i<dim; i++){
		for (int j(0); j<i; j++){
			cout<<0.0<<" " ;
		}
		cout << D[i] <<" " ;
		for (int j(i+1); j<dim; j++){
			cout<<0.0<<" " ;
		}
		cout <<endl ;		
	}	
	cout <<endl ;
	cout<<endl ;
	cout << "Solution:\n-----------"<<endl ;
	for (int i(0); i<dim; i++){
		cout << x[i] <<" " ;
	}
	cout <<endl ;
	cout<<endl ;
}

void Factor::solve(){
	double r ;
	// Solve Lx=b
	for (int i(0); i<dim; i++){
		r = b[i] ;
		for (int j(0); j<i; j++){
			r -= L[i][j]*x[j] ;
		}
		x[i] = r ;
	}

	// solve Dx=x
	for (int i(0); i<dim; i++){
		x[i] /= D[i] ;
	}

	// solve Lt x = x
	for (int i(dim-1); i>=0; i--){
		r = x[i] ;
		for (int j(i+1); j<dim; j++){
			r -= L[j][i]*x[j] ;
		}
		x[i] = r ;
	}
}