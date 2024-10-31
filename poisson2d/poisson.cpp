#include "tools.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
using namespace std ;

#define pi M_PI

class Poisson{
	private:
		int N ;
		int M ;
		double ax ;
		double by ;
		double dx ;
		double dy ;
		double** S ;
		double** U ;
		double* d ;
	public:
		Poisson() ;
		~Poisson() ;
		void start() ;
		void solve() ;
		void finalize() ;	
} ;

int main(){
	char choice = 'y' ;
	while (choice == 'y'){
		system("cls") ;
		Poisson p ;
		p.start() ;
		p.solve() ;
		p.finalize() ;
		system("python poisson.py") ;
		cout << "Reset ? [y|n]" ; cin >> choice ; 	
	}

	return 0 ;
}

Poisson::Poisson(){
		Gen g ;
		g.generate() ;
		N = g.getN() ;
		M = g.getM() ;
		ax = g.getax() ;
		by = g.getby() ;
		dx = g.getdx() ;
		dy = g.getdy() ;
		double a = -pow(dy/dx,2) ;
		try{
			double** F = g.getF() ;
													
			S = new double*[N] ;
			for (int i=0; i<N; i++){
				S[i] = new double[M] ;
			}
													
			// ligne 1
			S[0][0] = F[1][1]+F[0][1]-a*F[1][0] ;
			for (int j=1; j<M-1; j++){
				S[0][j] = F[1][j+1]+F[0][j+1] ;
			}
			S[0][M-1] = F[1][M]+F[0][M]-a*F[1][M+1] ;
													
			//ligne 2 à N-1
			for (int i=1; i<=N-2; i++){
				S[i][0] = F[i+1][1]-a*F[i][0] ;
				for (int j=1; j<M-1; j++){
					S[i][j] = F[i+1][j+1] ;
				}
				S[i][M-1] = F[i+1][M]-a*F[i+1][M+1] ;
			}
													
			//ligne N
			S[N-1][0] = F[N][1]-a*F[N][0]+F[N+1][1] ;
			for (int j=1; j<M-1; j++){
				S[N-1][j] = F[N][j+1]+F[N+1][j+1] ;
			}
			S[N-1][M-1] = F[N][M]-a*F[N][M+1]+F[N+1][M] ;
													
			// valeurs propres
			d = new double[M] ;
			for (int j=1; j<=M; j++){
				d[j-1] = 2*(1-a*(1-cos(pi*j/(M+1)))) ;
			}
													
			U = new double*[N] ;
			for (int i=0; i<N; i++){
				U[i] = new double[M] ;
			}				
		}
		catch(exception& e){				
			cout << "Erreur mémoire" << e.what() << endl ;
		}
}

Poisson::~Poisson(){
	for (int i=0; i<M; i++){
		delete[] S[i] ;
	}
	for (int i=0; i<N; i++){
		delete[] U[i] ;
	}
	delete[] S ;
	delete[] U ;
	delete[] d ;
}

void Poisson::start(){
	for (int j=0; j<N; j++){
		DST(S[j], M) ;
		for (int k=0; k<M; k++){
			S[j][k] *= sqrt(2.0/(M+1)) ;
		}
	}
}

void Poisson::solve(){
	try{
		double** T = new double*[3] ;
		for (int i=0; i<3; i++){
			T[i] = new double[N] ;
		}
		T[0][0] = 0.0 ;
		for (int j=1; j<N; j++){
			T[0][j] = -1 ;
			T[2][j-1] = -1 ;
		}
		T[2][N-1] = 0.0 ;
									
		double* u = new double[N] ;
		double* h = new double[N] ;
		for (int j=0; j<M; j++){
			for (int k=0; k<N; k++){
				h[k] = S[k][j] ;
				T[1][k] = d[j] ;
			}
			tridiag(T,h,u,N) ;
			for (int k=0; k<N; k++){
				U[k][j] = u[k] ;
			}
		}
		for (int i=0; i<3; i++){
			delete[] T[i] ;
		}
		delete[] T ;
		delete[] u ;
		delete[] h ;
	}
	catch(exception& e){
		cout << "Erreur memoire"<<e.what()<<endl ;
	}
}

void Poisson::finalize(){
	for (int i=0; i<N; i++){
		DST(U[i],M) ;
	}
	ofstream f("solution.txt") ;
	double x, y ;
	for (int i=0; i<N; i++){
		for (int j=0; j<M; j++){
			x = ax+(dx*j) ;
			y = by-(dy*i) ;
			f << x <<" "<< y <<" " << abs(U[i][j]-u(x,y)) << endl ;
		}
	}
	f.close() ;
}