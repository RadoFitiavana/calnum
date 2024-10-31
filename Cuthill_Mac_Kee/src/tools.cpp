#include <vector>
#include "tools.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <graphviz/gvc.h>

using namespace std ;

/************************************** Usefull functions **************************************/
vector<int> is_in(int c, vector<int>& v){
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
		if (is_in(c,container)[0] == 0){
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
	vector<int> tmp = is_in(c, container) ;
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
/***************************************************************************************************************/

/********************************************* The Cuthill Mac Kee processes *********************************/

Cuthill::Cuthill(){
	ifstream f("matrix.txt") ;
	if (f){
		f >> n ;
		float val ;
		vector<int> tmp ;
		for (int i(0); i<n; i++){
			for (int j(0); j<i; j++){
				f >> val ;
				if (val != 0){
					tmp.push_back(j) ;
				}
			}
			 f >> val ;
			for (int j(i+1); j<n; j++){
				f >> val ;
				if (val != 0){
					tmp.push_back(j) ;
				}
			}
			L.push_back(tmp) ;
			tmp.clear() ;
		}
	}

	else{
		cerr << "Can't open file matrix.txt" << endl ;
	}
}


vector<int> Cuthill::nei(vector<int>& E){
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

vector<int> Cuthill::ex(int s){
	int e(0) ;
	vector<int> l_e ;
	vector<int> a = {s} ;
	vector<int> b = L[s] ;
	vector<int> c ;
	vector<int> tmp ;
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
	l_e.push_back(a[0]) ;
	l_e.push_back(a.size()) ;
	return l_e ;
}

int Cuthill::get_size(){
	return n ;
}

vector<int> Cuthill::get_nei(int s){
	return L[s] ;
}

int Cuthill::start_node(int s){
	int N(this->n) ;
	vector<int> prev ;
	int start(s) ;
	vector<int> cur = this->ex(start);
	N -= 1 ;
	while (N > 0) {
		prev = cur ;
		cur = this->ex(prev[1]) ;
		N -= prev[2] ;
		if (cur[0] > prev[0]){
			start = cur[1] ;
		}
	}
	cout << "start node: "<< start <<endl ;
	return start ;
}

vector<int> Cuthill::cmk(int start){
	vector<int> nodes = {start} ;
	vector<int> to_visit = L[start] ;
	int next ;
	int i(0) ;
	int pos ;
	while (nodes.size() < this->n){
		// sort the list of to_visit
		while (to_visit.size() > 0){
			next = to_visit[0] ;
			pos = 0 ;
			for (int i(0); i<to_visit.size(); i++){
				L[to_visit[i]] = L[to_visit[i]] - nodes ;
				if (L[to_visit[i]].size() < L[next].size()){
					next = to_visit[i] ;
					pos = i ;				
				}
			}
			nodes.push_back(next) ;
			to_visit.erase(to_visit.begin()+pos) ;
		}
		i += 1 ;
		to_visit = L[nodes[i]] - nodes ;
	}
	reverse(nodes.begin(), nodes.end());
	return nodes ;
}

/************************************************************************************************/

/*******************************Range the profil of the matrix***************************************/

void range(double** A, std::vector<double>& b, std::vector<int>& perm) {
    int n = perm.size();
    std::vector<std::vector<double>> temp(n, std::vector<double>(n, 0.0));
    std::vector<double> temp_b(n);

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            temp[i][j] = A[perm[i]][perm[j]];
            temp[j][i] = temp[i][j];
        }
        temp_b[i] = b[perm[i]];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = temp[i][j];
        }
        b[i] = temp_b[i];
    }
}


Solver::Solver(){
	int start ;
	cout << "Initialization for Cuthill: "; cin >> start ;
	cout << endl ;
	Cuthill C ;
	perm = C.cmk(C.start_node(start)) ;
	ifstream f("matrix.txt") ;
	if (f){
		int dim ;
		f >> dim ;
		double** tmp(nullptr) ;
		try{
			tmp = new double*[dim] ;
			for (int i(0); i<dim; i++){
				tmp[i] = new double[dim] ;
				for (int j(0); j<dim; j++){
					f >> tmp[i][j] ;
				}
			}
			double s ;
			for (int i(0); i<dim; i++){
				f >> s ;
				b.push_back(s) ;
			}

			Graph g1(tmp,dim) ;
			g1.render_graph("dot -Tpng graph.dot -o initial.png") ;

			range(tmp,b,perm) ;

			Graph g2(tmp,dim) ;
			g2.render_graph("dot -Tpng graph.dot -o arranged.png") ;

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

Solver::~Solver(){
	delete[] LP ;
}

void Solver::display(){
	cout << "A'P:\n----------"<<endl ;
	for (int i(0); i<AP.size(); i++){
		cout << AP[i] <<" " ;
	}
	cout<<endl ;
	cout<<endl ;	
	cout << "L'P:\n----------"<<endl ;
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

void Solver::factorize(){
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

void Solver::solve(){
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
	 vector<double> tmp(dim) ;
	for (int i(0); i<x.size(); i++){
		tmp[perm[i]] = x[i] ; 
	}
	for (int i(0); i<dim; i++){
		x[i] = tmp[i] ;
	}	
}

Graph::Graph(double** A, int size){
    for (int i(0); i<size; i++){
        for (int j(0); j<i; j++){
            if (A[i][j] != 0){
                this->edges.push_back({i,j}) ;
            }
        }
    }
}

void Graph::render_graph(const char* command){
    ofstream f("graph.dot") ;
    if (f){
        f << "graph G{"<<endl ;
        for (int i(0); i<this->edges.size(); i++){
            f << "\t" << this->edges[i][0]<<" -- "<<this->edges[i][1]<< " ;" <<endl ;
        }
        f << "}" ;
        f.close() ;
        system(command) ;
        system("rm graph.dot") ;
    }
    else{
        cerr<< "Can't create graph" << endl ;
    }
}

