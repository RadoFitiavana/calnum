#ifndef RANGE_H
#define RANGE_H

#include <vector>

class Profil {
	private:
		std::vector<double> AP ; // The profil of the matrix A
		double* LP ; // the profil of L from the LDLt decomposition of A
		std::vector<int> nDiag ;  // the vector of the position of each diagonal into the profil
		std::vector<int> p ; // the vector of the indicies of the firsts non zero element of lines of matrix 
		std::vector<double> b ; // the second member of the linear equation to solve
		std::vector<double> x ; // the container of the final solution of the system
		std::vector<double> diag ; // the diagonal of the diagonal matrix D from the LDLt decomposition of A

	public:
		Profil() ; // constructor
		~Profil() ; // destructor
		void factorize() ; // LDLt factorization 
		void solve() ; // solve the equation Ax=b
		void display() ; // display result (LDLt + solution)
} ;

#endif