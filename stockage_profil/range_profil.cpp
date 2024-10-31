#include "range_profil.h"
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std ;

Profil::Profil(){
	ifstream f("matrix.txt") ;
	if (f){
		int dim ;
		f >> dim ;
		double** tmp ;
		try{
			tmp = new double*[dim] ;
			for (int i(0); i<dim; i++){
				tmp[i] = new double[i] ;
				for (int j(0); j<=i; j++){
					f >> tmp[i][j] ;
				}
			}

			nDiag.push_back(0) ;
			AP.push_back(tmp[0][0]) ;
			p.push_back(0) ;
			int pos(1) ;
			int j ;
			for (int i(1); i<dim; i++){
				j = 0 ;
				while (j<i){
					if (tmp[i][j]!=0){
						break ;
					}
					j += 1 ;
				}
				if (j == i-1){
					p.push_back(i) ;
				}
				else{
					p.push_back(j) ;
					for (int k(j); k<i; k++){
						AP.push_back(tmp[i][k]) ;
						pos += 1 ;
					}
				}
				AP.push_back(tmp[i][i]) ;
				nDiag.push_back(pos) ;
				pos += 1 ;
			}

			for (int i(0); i<dim; i++){
				delete[] tmp[i] ;
			}
			delete[] tmp ;

			LP = new double[AP.size()] ;

			double s ;
			for (int i(0); i<dim; i++){
				f >> s ;
				b.push_back(s) ;
			}

		}
		catch(exception& e){
			cerr << "Memory error: " <<e.what() << endl ;
		}
		f.close() ;
	}
	else{
		cerr << "Can't open file matrix.txt" <<endl ;
	}
}

Profil::~Profil(){
	delete[] LP ;
}

void Profil::display(){
	cout << "AP:\n----------"<<endl ;
	for (int i(0); i<AP.size(); i++){
		cout << AP[i] <<" " ;
	}
	cout<<endl ;
	cout<<endl ;	
	cout << "LP:\n----------"<<endl ;
	for (int i(0); i<AP.size(); i++){
		cout << LP[i] <<" " ;
	}
	cout<<endl ;
	cout<<endl ;
	cout << "nDiag:\n----------"<<endl ;
	for (int i(0); i<nDiag.size(); i++){
		cout << nDiag[i]+1 << " " ;
	}
	cout<<endl ;
	cout <<endl ;
	cout << "pi:\n----------"<<endl ;
	for (int i(0); i<p.size(); i++){
		cout << p[i]+1 << " " ;
	}
	cout<<endl ;
	cout <<endl ;
	cout << "diag:\n----------"<<endl ;
	for (int i(0); i<diag.size(); i++){
		cout << diag[i] <<" " ;
	}
	cout<<endl ;
	cout<<endl ;
	cout <<"solution:\n----------"<<endl ;
	for (int i(0); i<x.size(); i++){
		cout << x[i] <<" " ;
	}
	cout<<endl ;
	cout<<endl ;	
}

void Profil::factorize(){
	int dim = nDiag.size() ;
	diag.push_back(AP[0]) ;
	LP[0] = 1 ;
	double r ;
	int max ;
	for (int i(1); i<dim; i++){
		LP[nDiag[i]] = 1 ;

		for (int j(p[i]); j<i; j++){
			r = 0.0 ;
			if (p[i]<=p[j]){
				max = p[j] ;
			}
			else{
				max = p[i] ;
			}
			for (int k(max); k<j; k++){
				r += LP[nDiag[i]-i+k]*diag[k]*LP[nDiag[j]-j+k] ;
			}
			LP[nDiag[i]-i+j] = (AP[nDiag[i]-i+j] - r)/diag[j] ;
		}

		r = 0.0 ;
		for (int k(p[i]); k<i; k++){
			r += diag[k]*pow(LP[nDiag[i]-i+k],2) ;
		}

		diag.push_back(AP[nDiag[i]] - r) ;
	}
}

void Profil::solve(){
	int dim = nDiag.size() ;
	double r ; 
	// solve Lx=b
	for (int i(0); i<dim; i++){
		r = b[i] ;
		for (int j(p[i]); j<i; j++){
			r -= LP[nDiag[i]-i+j]*x[j] ;
		}
		x.push_back(r) ;
	}

	// solve Dx=x
	for (int i(0); i<dim; i++){
		x[i] /= diag[i] ;
	}

	// solve Lt x = x
	for (int i(dim-1); i>=0; i--){
		r = x[i] ;
		for (int j(i+1); j<dim; j++){
			if (i >= p[j]){
				r -= LP[nDiag[j]-j+i]*x[j] ;
			}
		}
		x[i] = r ;
	}	
}