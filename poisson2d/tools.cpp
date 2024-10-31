#include "tools.h"
#include <complex>
#include <iostream>
#include <cmath>
using namespace std ;

#define pi M_PI

void range(double* d, int n){
	double* t = new double[n] ;
	int p = n/2 ;
	for (int i=0; i<p; i++){
		t[i] = d[2*i] ;
		t[i+p] = d[2*i + 1] ;
	}
	for (int i=0; i<n; i++){
		d[i] = t[i] ;
	}
	delete[] t ;
}
void fft(double* d, complex<double>* z, int n){
	complex<double> i(0.0, 1.0) ;
	if (n == 2){
		z[0] = d[0] + d[1] ;
		z[1] = d[0] - d[1] ;
	}
	else{
		int p = n/2 ;
		range(d, n) ;
		fft(d, z, p) ;
		fft(d+p, z+p, p) ;
		complex<double> a ;
		for (int k=0; k<p; k++){
			z[k+p] *= exp(i*(pi*k/p)) ;
			a = z[k] ;
			z[k] += z[k+p] ;
			z[k+p] = a - z[k+p] ;
		}
	}
}

void fst(double* v, int n){
	double* d = new double[n] ;
	complex<double>* z = new complex<double>[n] ;
	d[0] = 0.0 ;
	for (int k=1; k<n; k++){
		d[k] = 0.5*(v[k]-v[n-k])+sin(k*pi/n)*(v[k]+v[n-k]) ;
	}
	fft(d, z, n) ;
	delete[] d ;
	for (int j=0; j<n; j++){
		z[j] *= 2.0/n ;
	}
	double r = sqrt(n/2.0) ;
	v[1] = r*(z[0].real()) ;
	for (int j=1; j<n/2; j++){
		v[2*j] = 2*r*(z[j].imag()) ;
		v[2*j+1] = v[2*j-1] + 2*r*(z[j].real()) ;
	}
	delete[] z ;
}

void DST(double* v, int n){
	int p = n+1 ;
	double* Sv = new double[p] ;
	Sv[0] = 0.0 ;
	for (int j=1; j<p; j++){
		Sv[j] = v[j-1] ;
	}
	fst(Sv, p) ;
	for (int j=1; j<p; j++){
		v[j-1] = Sv[j] ;
	}
}

void tridiag(double** A, double* b, double* x, int n){
	double k ;
	for (int j=0; j<n-1; j++){
		k = A[2][j]/A[1][j] ;
			A[1][j+1] -= k*A[0][j+1] ;
			b[j+1] -= k*b[j] ;
	}
	x[n-1] = b[n-1]/A[1][n-1] ;
	for (int j=n-2; j>=0; j--){
		x[j] = (b[j] - A[0][j+1]*x[j+1])/A[1][j] ;
	}
}

Gen::Gen(){
	cout << "N = " ; cin >> N ;
	cout << "M = " ; cin >> M ;
	cout << "ax = " ; cin >> ax ;
	cout << "bx = " ; cin >> bx ;
	cout << "ay = " ; cin >> ay ;
	cout << "by = " ; cin >> by ;
	try{
		F = new double*[N+2] ;
		for (int i=0; i<N+2; i++){
			F[i] = new double[M+2] ;
		}
	}
	catch(exception& e){
		cout << "Erreur memoire" << e.what() << endl ;
	}
}

Gen::~Gen(){
	for (int i=0; i<N+2; i++){
		delete[] F[i] ;
	}
	delete[] F ;
}

void Gen::generate(){
	double dx = (bx-ax)/(M+1) ;
	double dy = (by-ay)/(N+1) ;
	F[0][0] = g(ax,by) ;
	F[0][M+1] = g(bx,by) ;
	F[N+1][0] = g(ax,ay) ;
	F[N+1][M+1] = g(bx,ay) ;
	for (int j=1; j<=M; j++){
		F[0][j] = g(ax+dx*j,by) ;
		F[N+1][j] = g(ax+dx*j, ay) ;
	}
	for (int i=1; i<=N; i++){
		F[i][0] = g(ax, by-dy*i) ;
		F[i][M+1] = g(bx, by-dy*i) ;
	}
	for (int i=1; i<=N; i++){
		for (int j=1; j<=M; j++){
		F[i][j] = -pow(dy,2)*f(ax+(dx*j),by-(dy*i)) ;
		}
	}
}

int Gen::getN(){
	return N ;
}

int Gen::getM(){
	return M ;
}

double Gen::getdx(){
	return (bx - ax)/(M+1) ;
}

double Gen::getdy(){
	return (by - ay)/(N+1) ;
}

double Gen::getax(){
	return ax ;
}

double Gen::getbx(){
	return bx ;
}

double Gen::getay(){
	return ay ;
}

double Gen::getby(){
	return by ;
}

double** Gen::getF(){
	return F ;
}

double Gen::g(double x, double y){
	//double r = pow(x,2) + pow(y,2) ;
	//return exp(-0.5*r) ;
	return x*x + y*y ;
}

double Gen::f(double x, double y){
	//return a*g(x,y) ;
	return 4 ;
}

double u(double x, double y){
	/*double r = pow(x,2) + pow(y,2) ;
	return exp(-0.5*r) ;*/
	return x*x + y*y ;	
}