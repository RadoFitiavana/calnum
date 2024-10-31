#ifndef POISSON_H
#define POISSON_H

#include <complex>
#include <cmath>
#include <omp.h>
#include <thread>
#include <stdexcept>
#include <functional>
#include <type_traits>
#include <vector>
#include <fstream>
#include <string>
using namespace std ;

namespace Poisson {

static int const nCores = std::thread::hardware_concurrency() ;

class ParallSetup {
	private :
		static int nThreads ;
		static bool PARALLEL ;

	public :
		static void setnTreads(int n) {
			nThreads = ((n < nCores)? n : nCores) ; 
		}
		static int getnThreads() {
			return nThreads ;
		}
		static void setPARALLEL(bool const enable) {
	    	PARALLEL = enable ;
    		if (PARALLEL) {
				omp_set_num_threads(nThreads) ;
        		omp_set_dynamic(1); // Enable dynamic adjustment of the number of thread
    		}	 
			else {
        		omp_set_num_threads(1); // Disable parallel computing by setting to one thread
        		omp_set_dynamic(0); // Disable dynamic adjustment
    		}		
		}
		static int getPARALLEL() {
			return PARALLEL ;
		}
} ;
bool ParallSetup::PARALLEL = false ;
int ParallSetup::nThreads = nCores/2 ; 

template <typename T>
class FDST {
	static_assert(
		is_same<T,float>::value || is_same<T,double>::value,
		"Allowed number type is either float or double \n"
	) ;
    private :
        int N ;
        T *signal ;
        void swapEvenOdd(T*, int) ;
        void fft(T*, std::complex<T>*, int) ;
        static T const pi ;

    public :
        FDST() ;
        FDST(T*, int) ;
        void setSignal(T*, int) ;
        void runFDST() ;

} ;
template <typename T>
T const FDST<T>::pi = (T)M_PI ;

template <typename T>
FDST<T>::FDST() {signal(nullptr) ;}

template <typename T>
FDST<T>::FDST(T *signal, int N) {
    this->signal = signal ;
    this->N = N ;
}

template <typename T>
void FDST<T>::swapEvenOdd(T* V, int P) {
	T* t = new T[P] ;

	if (ParallSetup::getPARALLEL()) {
		////omp_set_num_threads
		#pragma omp parallel
		{
			#pragma omp for
			for (int i=0; i<P/2; i++){
				t[i] = V[2*i] ;
				t[i+P/2] = V[2*i + 1] ;	
			}
		}
		//omp_set_num_threads(((P < nCores)? P : nCores)) ;
		#pragma omp parallel
		{
			#pragma omp for
			for (int i=0; i<P; i++){
				V[i] = t[i] ;
			}
		}
	}
	else {
		for (int i=0; i<P/2; i++){
			t[i] = V[2*i] ;
			t[i+P/2] = V[2*i + 1] ;	
		}
		for (int i=0; i<P; i++){
			V[i] = t[i] ;
		}		
	}

	delete[] t ;
}

template <typename T>
void FDST<T>::fft(T *V, complex<T> *fV, int P) {
	complex<T> i(0.0, 1.0) ;
	if (P == 2){
		fV[0] = V[0] + V[1] ;
		fV[1] = V[0] - V[1] ;
	}
	else{
		swapEvenOdd(V, P) ;

		if(ParallSetup::getPARALLEL()) {
			//omp_set_num_threads(2) ;
			#pragma omp parallel
			{
				#pragma omp single
				{
					#pragma omp task
					{
						fft(V, fV, P/2) ;
					}
					#pragma omp task 
					{
						fft(V+P/2, fV+P/2, P/2) ;
					}

					#pragma omp taskwait
				}
			}
		}
		else {
			fft(V, fV, P/2) ;
			fft(V+P/2, fV+P/2, P/2) ;
		}

		complex<T> a ;

		for (int k=0; k<P/2; k++){
			fV[k+P/2] *= exp(i*(2*k*pi/P)) ;
			a = fV[k] ;
			fV[k] += fV[k+P/2] ;
			fV[k+P/2] = a - fV[k+P/2] ;
		}
	}
}

template <typename T>
void FDST<T>::runFDST() {
    int n = N+1 ;
	T* d = new T[n] ;
	complex<T>* z = new complex<T>[n] ;
	d[0] = 0.0 ;

	if(ParallSetup::getPARALLEL()) {
		//omp_set_num_threads(((n < nCores+1)? n-1 : nCores)) ;
		#pragma omp parallel
		{
			#pragma omp for
			for (int k=1; k<n; k++){
				d[k] = 0.5*(signal[k-1]-signal[N-k])+sin(k*pi/n)*(signal[k-1]+signal[N-k]) ;
			}
		}
	}
	else {
		for (int k=1; k<n; k++){
			d[k] = 0.5*(signal[k-1]-signal[N-k])+sin(k*pi/n)*(signal[k-1]+signal[N-k]) ;
		}
	}

	fft(d, z, n) ;
	delete[] d ;
	d = nullptr ;

	if(ParallSetup::getPARALLEL()) {
		//omp_set_num_threads(((n < nCores)? n : nCores)) ;
		#pragma omp parallel
		{
			#pragma omp for
			for (int j=0; j<n; j++){
				z[j] *= 2.0/n ;
			}
		}
	}
	else {
		for (int j=0; j<n; j++){
			z[j] *= 2.0/n ;
		}		
	}
	
	T r = sqrt(n/2.0) ;
	signal[0] = r*(z[0].real()) ;
	for (int j=1; j<n/2; j++){
		signal[2*j-1] = 2*r*(z[j].imag()) ;
		signal[2*j] = signal[2*j-2] + 2*r*(z[j].real()) ;
	}
	delete[] z ;
}

template <typename T>
void tridiagSolver(int N, T **A, T *b, T *x) {
	if(ParallSetup::getPARALLEL()) {
		//omp_set_num_threads(((N < nCores+1)? N-1 : nCores)) ;
		#pragma omp parallel
		{
			#pragma omp for
			for (int j=0; j<N-1; j++){
				A[1][j+1] -= (A[2][j]/A[1][j])*A[0][j+1] ;
				b[j+1] -= (A[2][j]/A[1][j])*b[j] ;
			}	
		}

		x[N-1] = b[N-1]/A[1][N-1] ;

		#pragma omp parallel
		{
			#pragma omp for
			for (int j=N-2; j>=0; j--){
				x[j] = (b[j] - A[0][j+1]*x[j+1])/A[1][j] ;
			}	
		}
	}
	else {
		for (int j=0; j<N-1; j++){
			A[1][j+1] -= (A[2][j]/A[1][j])*A[0][j+1] ;
			b[j+1] -= (A[2][j]/A[1][j])*b[j] ;
		}
		x[N-1] = b[N-1]/A[1][N-1] ;
		for (int j=N-2; j>=0; j--){
			x[j] = (b[j] - A[0][j+1]*x[j+1])/A[1][j] ;
		}		
	}
}

template <typename T>
class Discretizer {
	static_assert(
		is_same<T,float>::value || is_same<T,double>::value,
		"Allowed number type is either float or double \n"
	) ;
	private :
		int N, M ;
		T ax, bx, ay, by ;
		function<T(T,T)> Lf ;
		function<T(T,T)> Bordf ;

	public :
		Discretizer() ;
		void setDomain(int, T) ;
		void setShape(int, int) ;
		vector<vector<T>> getData() ;
		vector<int> getShape() ;
		vector<T> getDomain() ;
		void setLaplacian(function<T(T,T)>) ;
		void setBorderFunction(function<T(T,T)>) ;
		T getLaplacian(T, T) ;
		T getBorderFunction(T, T) ;

} ;

template <typename T>
Discretizer<T>::Discretizer() {
	ax = -5.0 ; bx = 5.0 ;
	ay = -5.0 ; by = 5.0 ;
	Bordf = [](T x, T y) {
				return 0.0 ;
			} ;
	Lf = [](T x, T y) {
			T r = pow(x,2) + pow(y,2), a = pow(x,2)+pow(y,2)-2 ;
			return a*exp(-0.5*r) ;
		 } ;
	N = 15 ;
	M = 15 ;
}

template <typename T>
void Discretizer<T>::setDomain(int bound, T val) {
	switch(bound) {
		case 0 :
			ax = val ;
			break ;
		case 1 :
			bx = val ;
			break ;
		case 2 :
			ay = val ;
			break ;
		case 3 :
			by = val ;
			break ;
		default :
			throw runtime_error(
				"bound(The first parameter) must be one of the following: 0(ax),1(bx),2(ay),3(by)\n"
			) ;
			break ;
	}
}

template <typename T>
void Discretizer<T>::setShape(int axis, int size) {
	int k = (int) (log(size)/log(2)) ;
	bool status = (size+1 == (int)pow(2,k+1)) ;
	if (status) { 
		switch(axis) {
			case 0 :
				N = size ;
				break ;
			case 1 :
				M = size ;
				break ;
			default :
				throw runtime_error(
					"axis(The first parameter) must be one of the following: 0(x-axis), 1(y-axis)\n"
				) ;
				break ;
		}
	}
	else {
		throw runtime_error(
			"Allowed discretization stepsize is a positive integer of the form 2^k - 1\n"
		) ;
	}

}

template <typename T>
vector<int> Discretizer<T>::getShape() {
	return {N,M} ; 
}

template <typename T>
vector<T> Discretizer<T>::getDomain() {
	return {ax,bx,ay,by} ;
}

template <typename T>
void Discretizer<T>::setLaplacian(function<T(T,T)> Lf) {
	this->Lf = Lf ;
}

template <typename T>
T Discretizer<T>::getLaplacian(T x, T y) {
	return this->Lf(x,y) ;
}

template <typename T>
void Discretizer<T>::setBorderFunction(function<T(T,T)> Bordf) {
	this->Bordf = Bordf ;
}

template <typename T>
T Discretizer<T>::getBorderFunction(T x, T y) {
	return this->Bordf(x,y) ;
}

template <typename T>
vector<vector<T>> Discretizer<T>::getData() {
	T dx((bx-ax)/(M+1)), dy((by-ay)/(N+1)) ;
	vector<vector<T>> F(N+1,vector<T>(M+1)) ;
	F[0][0] = Bordf(ax,by) ;
	F[0][M+1] = Bordf(bx,by) ;
	F[N+1][0] = Bordf(ax,ay) ;
	F[N+1][M+1] = Bordf(bx,ay) ;

	if (ParallSetup::getPARALLEL()) {
		#pragma omp parallel
		{
			#pragma omp for
			{
				for (int j=1; j<=M; j++){
					F[0][j] = Bordf(ax+dx*j,by) ;
					F[N+1][j] = Bordf(ax+dx*j, ay) ;
				}
			}
		}
		#pragma omp parallel
		{
			#pragma omp for
			{
        		for (int i=1; i<=N; i++){
					F[i][0] = Bordf(ax, by-dy*i) ;
					F[i][M+1] = Bordf(bx, by-dy*i) ;
				}
			}
		}
		#pragma omp parallel
		{
			#pragma omp for
			{
				for (int i=1; i<=N; i++){
					for (int j=1; j<=M; j++){
						F[i][j] = -pow(dy,2)*Lf(ax+(dx*j),by-(dy*i)) ;
					}
				}
			}
		}
	}
	else {
		for (int j=1; j<=M; j++){
			F[0][j] = Bordf(ax+dx*j,by) ;
			F[N+1][j] = Bordf(ax+dx*j, ay) ;
		}
		for (int i=1; i<=N; i++){
			F[i][0] = Bordf(ax, by-dy*i) ;
			F[i][M+1] = Bordf(bx, by-dy*i) ;
		}
		for (int i=1; i<=N; i++){
			for (int j=1; j<=M; j++){
				F[i][j] = -pow(dy,2)*Lf(ax+(dx*j),by-(dy*i)) ;
			}
		}	
	}
	return F ;
}

template <typename T>
class Poisson2DSolver {
	static_assert(
		is_same<T,float>::value || is_same<T,double>::value,
		"Allowed number type is either float or double \n"
	) ;
	private :
		int N, M ;
		T **S, **U, *eigVal, ax, bx, ay, by ;
		void start() ;
		void finalize(bool) ;

	public :
		Poisson2DSolver() ;
		Poisson2DSolver(Discretizer<T>&) ;
		Poisson2DSolver(int, int, T, T, T, T, string) ;
		Poisson2DSolver(int, int, T, T, T, T, vector<vector<T>>&) ;
		void setDomain(int, T) ;
		void setShape(int, int) ;
		void runSolver(bool) ;
		void plot() ;		
		~Poisson2DSolver() ;
} ;

template <typename T>
Poisson2DSolver<T>::Poisson2DSolver() {
	S(nullptr); U(nullptr); eigVal(nullptr) ;
}

template <typename T>
Poisson2DSolver<T>::Poisson2DSolver(Discretizer<T>& dis) {
	
}

}

#endif