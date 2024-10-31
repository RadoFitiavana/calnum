#ifndef TOOLS_H
#define TOOLS_H

#include <vector>

// the function isin verify if the integer c is present
// into the vector v
std::vector<int> isin(int c, std::vector<int>& v) ;

void extends(std::vector<int>& a, std::vector<int>& b) ; // append the vector a to the vector b 

void append(int c, std::vector<int>& v) ; // append the integer c to the vector v

void del(int c, std::vector<int>& v) ; // delete the integer c into the vector v

// calcul C <- A-B and return C 
std::vector<int> operator-(std::vector<int>& A, std::vector<int>& B) ;

// range a given symetric matrix A according to a given permutation perm
void range(double** A, std::vector<double>& b, std::vector<int>& perm) ; 

class Cuthill {

private:
	int n ; // number of nodes of the graph
	std::vector<std::vector<int>> L; // List of the neigbors of each node

public:
	Cuthill() ; // constructor of the Cuthill
			  //the constructor read data from the file matrix.txt	

	// return the neigbors of the vector v
	std::vector<int> nei(std::vector<int>& v) ;

	// return a vector containning the excentricity of the node, the first node of the last stage and the size of the last stage 
	std::vector<int> ex(int s) ;

	// getter of the size of the graph
	int get_size() ;

	// return the neighbors of the node s
	std::vector<int> get_nei(int s) ;

	// return the start node for the Cuthill Mac Kee algorithm
	int start_node(int s) ;

	// return the configuration from Cuthill Mac Kee algorithm
	std::vector<int> cmk(int start) ;

} ;

class Solver {

	private:
		std::vector<int> perm ;
		std::vector<double> AP ; // The profil of the matrix A
		double* LP ; // the profil of L from the LDLt decomposition of A
		std::vector<int> nDiag ;  // the vector of the position of each diagonal into the profil
		std::vector<int> p ; // the vector of the indicies of the firsts non zero element of lines of matrix 
		std::vector<double> b ; // the second member of the linear equation to solve
		std::vector<double> x ; // the container of the final solution of the system
		std::vector<double> diag ; // the diagonal of the diagonal matrix D from the LDLt decomposition of A

	public:
		Solver() ;
		~Solver() ;
		void factorize() ;
		void solve() ;
		void display() ;
} ;

class Graph{
    private:
        std::vector<std::vector<int>> edges ; 

    public:
        Graph(double**, int) ;
        void render_graph(const char*) ;

} ;

#endif