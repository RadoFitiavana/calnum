#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>

using namespace std ;

void span(int start, vector<int> to_visit, vector<int> path, vector<int> parents, vector<vector<int>>& nei, vector<vector<int>>& paths){

	if (to_visit.size() == 0){
		paths.push_back(path) ;
	}
	else{

		for (int i(0); i<to_visit.size(); i++){

			vector<int> new_parents = parents ;
			new_parents.push_back(start) ;

			int new_start(to_visit[i]) ;

			vector<int> new_path = path ;
			new_path.push_back(to_visit[i]) ;

			vector<int> new_to_visit = nei[new_start] ;

			for (int j(0); j<new_parents.size(); j++){
				for (int i(0); i<new_to_visit.size(); i++){
					if (new_to_visit[i] == new_parents[j]){
						new_to_visit.erase(new_to_visit.begin()+i) ;
						break ;
					}
				}
			}
			span(new_start,new_to_visit,new_path,new_parents,nei,paths) ;
		}
	}
}

vector<vector<int>> Span(int start, vector<vector<int>>& nei){
	vector<vector<int>> paths ;
	span(start,nei[start],{start},{},nei,paths) ;
	return paths ;
}

void depth(int start, vector<int> to_visit, vector<int> path, vector<int> prev, vector<vector<int>>& nei, vector<vector<int>>& paths){

	if (to_visit.size() == 0){
		paths.push_back(path) ;
	}
	else{

		for (int i(0); i<to_visit.size(); i++){

			vector<int> new_prev = prev ;
			new_prev.push_back(start) ;
			for (int k(0); k<nei[start].size(); k++){
				new_prev.push_back(nei[start][i]) ;
			}

			int new_start(to_visit[i]) ;

			vector<int> new_path = path ;
			new_path.push_back(to_visit[i]) ;

			vector<int> new_to_visit = nei[new_start] ;

			for (int j(0); j<new_prev.size(); j++){
				for (int i(0); i<new_to_visit.size(); i++){
					if (new_to_visit[i] == new_prev[j]){
						new_to_visit.erase(new_to_visit.begin()+i) ;
						break ;
					}
				}
			}
			depth(new_start,new_to_visit,new_path,new_prev,nei,paths) ;
		}
	}
}

vector<vector<int>> Depth(int start, vector<vector<int>>& nei){
	vector<vector<int>> paths ;
	depth(start,nei[start],{start},{},nei,paths) ;
	return paths ;
}

vector<int> intersect(vector<int>& u, vector<int> v){
	vector<int> u_v ;
	for (int i(0); i<u.size(); i++){
		for (int j(0); j<v.size(); j++){
			if (v[j] == u[i]){
				u_v.push_back(u[i]) ;
				v.erase(v.begin()+j) ;
				break ;
			}
		}
	}
	return u_v ;
}

vector<vector<int>> Connex(vector<vector<int>>& nei){
	vector<vector<int>> Paths ;
	vector<vector<int>> tmp = Depth(0,nei) ;
	for (int i(0); i<tmp.size(); i++){
		if (tmp[i].size() == nei.size()){
			Paths.push_back(tmp[i]) ;
		}
	}
	tmp = Paths ;
	Paths = {} ;
	vector<int> swap ;
	for (int i(0); i<tmp.size(); i++){
		for (int j(i+1); j<tmp.size(); j++){
			if (intersect(tmp[i],tmp[j]).size() == 1){
				for (int k(tmp[j].size()-1); k>0; k--){
					swap.push_back(tmp[j][k]) ;
				}
				for (int k(0); k<tmp[i].size(); k++){
					swap.push_back(tmp[i][k]) ;
				}
				Paths.push_back(swap) ;
				swap = {} ;
			}
		}
	}
	return Paths ;
}

void read_graph(vector<vector<float>>& mat, vector<vector<int>>& nei){
	ifstream f("graph.txt") ;
	if (f){
		int n ;
		f >> n ;
		vector<float> in_mat ;
		vector<int> in_nei ;
		float val ;
		for (int i(0); i<n; i++){
			for (int j(0); j<n; j++){
				f >> val ;
				in_mat.push_back(val) ;
				if (val != 0){
					in_nei.push_back(j) ;
				}
			}
			mat.push_back(in_mat) ;
			nei.push_back(in_nei) ;
			in_mat = {} ;
			in_nei = {} ;
		}
		f.close() ;
	}
	else{
		cerr << "Failed openning graph.txt" <<endl ;
	}
}

float path_weight(vector<int> path, vector<vector<float>> mat){
	float len(0) ;
	for (int i(0); i<path.size()-1; i++){
		len += mat[path[i]][path[i+1]] ;
	}
	return len ;
}

int main(){
	vector<vector<float>> mat ;
	vector<vector<int>> nei ;
	read_graph(mat,nei) ;
	int n(mat.size()) ;
	vector<vector<int>> all_paths ;
	//system("cls") ;
	/*vector<vector<int>> paths = Connex(nei) ;
	cout << paths.size() <<endl ;
	for (int i(0); i<paths.size(); i++){
		for (int j(0); j<paths[i].size(); j++){
			cout << paths[i][j] << " " ;
		}
		cout << endl ;
	}*/
	for (int start(0); start<n; start++){
		vector<vector<int>> paths = Depth(start,nei) ;
		
		cout << paths.size()<<endl ;
		for (int i(0); i<paths.size(); i++){
			if (paths[i].size() == n){
				all_paths.push_back(paths[i]) ;
				/*for (int j(0); j<paths[i].size(); j++){
					cout << paths[i][j] << " " ;
				}
				cout << endl ;*/
			}
		}
	}
	vector<int> shortest ;
	float min = path_weight(all_paths[0],mat) ;
	float cur ;
	for (int i(1); i<all_paths.size(); i++){
		cur = path_weight(all_paths[i],mat) ;
		if (cur < min){
			shortest = all_paths[i] ;
			min = cur ;
		}
	}
	float total(0) ;
	for (int i(0); i<mat.size(); i++){
		for (int j(0); j<i; j++){
			total += mat[i][j] ;
		}
	}
	cout << "total cost: " ;
	cout << total << endl ;
	cout << "Max saving: " << (total-min) << endl ;
	for (int i(0); i<shortest.size(); i++){
		cout << shortest[i] << " " ;
	}
	cout <<endl ;
	cout <<endl ;
	return 0 ;
}
