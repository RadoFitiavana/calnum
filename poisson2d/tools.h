#ifndef TOOLS_H
#define TOOLS_H

#include <complex>
using namespace std ;

void range(double*, int) ;
void fft(double*, complex<double>*, int) ;
void fst(double*, int) ;
void DST(double*, int) ;
void tridiag(double**, double*, double*, int) ;
double u(double, double) ;

class Gen{
private:
	int N ;
	int M ;
	double ax ;
	double bx ;
	double ay ;
	double by ;
	double** F ;
public:
	Gen() ;
	~Gen() ;
	void generate() ;
	int getN() ;
	int getM() ;
	double getdx() ;
	double getdy() ;
	double getax() ;
	double getay() ;
	double getbx() ;
	double getby() ;
	double** getF() ;
	double g(double, double) ;
	double f(double, double) ;
} ;

#endif