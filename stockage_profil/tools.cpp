#include <vector>
#include "tools.h"
#include <iostream>
#include <algorithm>
#include <fstream>
using namespace std ;

/********************Additional function**************************/
vector<int> isin(int c, vector<int>& v){
	int n = v.size() ;
	vector<int> res = {0,0} ;
	if ((c == v[0]) || (c == v[n-1])){
		res[0] = 1 ;
		if (c == v[0]){
			res[1] = 0 ;
		}
		else{
			res[1] = n-1 ;
		}
	}
	else{
		int min(0), max(n-1) ;
		int pos = (max+min)/2 ;
		while ((max-min) > 1){
			if (v[pos] == c){
				res[0] = 1 ;
				res[1] = pos ;
				break ;
			}
			if (v[pos] < c){
				min = pos ;
			}
			if (v[pos] > c){
				max = pos ;
			}
			pos = (max+min)/2 ;
		}
	}
	return res ;
}

void extends(vector<int>& v, vector<int>& set){
	int n = v.size() ;
	int c ;
	for (int i=0; i<n; i++){
		c = v[i] ;		
		append(c,set) ;
	}
}

void append(int c, vector<int>& container){
	if (container.size() != 0){
		if (isin(c,container)[0] == 0){
			int n = container.size()-1 ;
			if ((c < container[0]) || (c > container[n])){
				if (c < container[0]){
					container.insert(container.begin(),c) ;
				}
				else{
					container.push_back(c) ;
				}
			}
			else{
				int min(0), max(container.size()-1) ;
				int pos = (max+min)/2 ;
				while ((max-min)>1){
					if ((c<container[pos])&&(c>container[pos-1])){
						max = pos ;
						break ;
					}
					if ((c<container[pos+1])&&(c>container[pos])){
						max = pos+1 ;
						break ;
					}
					if (container[pos] < c){
						min = pos ;
					}
					if (container[pos] > c){
						max = pos ;
					}
					pos = (max+min)/2 ;
				}
				container.insert(container.begin()+max, c) ;
			}
		}	
	}
	else{
		container.push_back(c) ;
	}
}

void del(int c, vector<int>& container){
	int N(0) ;
	vector<int> tmp = isin(c, container) ;
	if (tmp[0] == 1){
		container.erase(container.begin()+tmp[1]) ;
	}
}

vector<int> operator-(vector<int>& A, vector<int>& B){
	vector<int> v = A ;
	int B_n = B.size() ;
	int c ;
 	for (int i=0; i<B_n; i++){
 		c = B[i] ;
 		del(c,v) ;
 	}
 	return v ;
}

void show(vector<int>& container){
	int n = container.size() ;
	cout <<"{" ;
	for (int i=0; i<n-1; i++){
		cout << container[i]+1 <<"," ;
	}
	if(n!=0){
		cout <<container[n-1]+1 ;
	}
	cout << "}" ;
}
/****************************************************************/

/*********************Graph methods definitions**********************/
Graph::Graph(){
	ifstream f("graph.txt") ;
	if (f){
		int N(0) ;
		int n_arrow ;
		vector<int> tmp(2) ;
		f >> n ;
		f >> n_arrow ;
		try{
			L = new vector<int>[n] ;
			for (int i(0); i<n_arrow; i++){
				for (int j(0); j<2; j++){
					f >> tmp[j] ;
					tmp[j] -= 1 ;
				}
				A.push_back(tmp) ;
				append(tmp[0],L[tmp[1]]) ;
				append(tmp[1],L[tmp[0]]) ;
			}
		}
		catch(exception& e){
			cerr << "An error occured while allocating: "<< e.what() << endl ;
		}
	}
	else{
		cerr << "Can't open file graph.txt" << endl ;
	}
}

Graph::~Graph(){
	delete[] L ;
}

void Graph::display(){
	cout << "Nodes: { " ;
	for (int i(0); i<this->n; i++){
		cout << "("<<i+1<<") " ;
	}
	cout << "}" <<endl ;
	cout << "Arrows: [ " ;
	for (int i(0); i<A.size(); i++){
		show(A[i]) ;
		cout <<" " ;
	}
	cout << "]" <<endl ;
	cout << "Neighbors: [ " ;
	for (int i(0); i<n; i++){
		cout <<"N("<<i+1<<")=" ;
		show(L[i]) ;
		cout <<" " ;
	}
	cout <<"]"<<endl ;
}

vector<int> Graph::nei(vector<int>& E){
	int p = E.size() ;
	vector<int> N_E ;
	if (p > 0){
		for (int i(0); i<p; i++){
			if ((E[i]>=0) && (E[i]<n)){
				extends(L[E[i]],N_E) ;
			}
		}
		N_E = N_E - E ;
	}
	return N_E ;
}

vector<int> Graph::ex(){
	int e ;
	vector<int> l_e ;
	vector<int> a ;
	vector<int> b ;
	vector<int> c ;
	vector<int> tmp ;
	for (int k(0); k<n; k++){
		e = 0 ;
		a = {k} ;
		b = L[k] ;
		while(b.size()>0){
			e += 1 ;
			tmp = nei(b) ;
			c = tmp - a ;
			tmp.clear() ;
			a.clear() ;
			a = b ;
			b.clear() ;
			b = c ;
			c.clear() ;
		}			
		l_e.push_back(e) ;
		
	}
	return l_e ;
}

int Graph::dima(){
	return d ;
}

void Graph::display_ex(){
	vector<int> l = ex() ;
	for (int k(0); k<l.size(); k++){
		cout << "e("<<k+1<<")="<< l[k] <<endl ;
	}
} 
/********************************************************************/