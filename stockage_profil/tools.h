#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <string>

// the function isin verify if the integer c is present
// into the vector v
std::vector<int> isin(int c, std::vector<int>& v) ;

void show(std::vector<int>& v) ; // display the vector v

void extends(std::vector<int>& a, std::vector<int>& b) ; // append the vector a to the vector b 

void append(int c, std::vector<int>& v) ; // append the integer c to the vector v

void del(int c, std::vector<int>& v) ; // delete the integer c into the vector v

// calcul C <- A-B and return C 
std::vector<int> operator-(std::vector<int>& A, std::vector<int>& B) ;

class Graph {
	/* we define a customize graph data structure 
	for matrix profile configuration*/

private:
	std::vector<std::vector<int>> A ; // list of undirected arrows of the graph
	int n ; // number of nodes of the graph
	std::vector<int>* L; // List of the neigbors of each node

public:
	Graph() ; // constructor of the graph
			  //the constructor read data from the file graph.txt	

	~Graph() ; // destructor of the graph

	// return the neigbors of the vector v
	std::vector<int> nei(std::vector<int>& v) ;

	// return the excentricity of the node 
	//int ex(int node) ;

	// return excentricity for all nodes as a vector formats
	std::vector<int> ex() ;

	// return the diameter of the graph
	int diam() ;

	// display the graph
	void display() ;

	// display excentricity of all nodes
	void display_ex() ;

} ;

#endif