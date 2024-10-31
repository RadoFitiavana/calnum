#ifndef TOOLS_H
#define TOOLS_H

class Factor{
	private:
		int dim ;
		double** A ;
		double** L ;
		double* D ;
		double* b ;
		double* x ;
	public:
		Factor() ;
		~Factor() ;
		void factorize() ;
		void display() ;
		void solve() ;
}; 

#endif