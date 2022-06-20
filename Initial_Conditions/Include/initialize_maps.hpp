#ifndef __initialize_maps__
#define __initialize_maps__
#include <sstream>
#include <string>
#include <map>
#include <utility>

std::map<std::pair<int,int>,double> initialize_map_neigh_angles(){
	
	std::map<std::pair<int,int>,double> neigh_angles;
	
	std::ifstream infile("Infile/Cos_angles_nearest_neighbors.txt");		

	std::string line;
	while (std::getline(infile, line))
	{
	    std::istringstream iss(line);
	    double a, b; 
	    double c;
	    if (!(iss >> a >> b >> c)) { 
	    	break; 
		}	 	
	   	std::pair<int,int> pair{int(a),int(b)};
	   	neigh_angles[pair]=c;
	}

	return neigh_angles;
}

std::map<std::vector<int>,int> initialize_map_vec_neigh(){
	
	std::map<std::vector<int>,int> vec_neigh;
	
	std::ifstream infile("Infile/MC_moves.txt");		

	std::string line;
	while (std::getline(infile, line))
	{
	    std::istringstream iss(line);
	    int a, b, c, d; 
	    if (!(iss >> a >> b >> c >>d)) { 
	    	break; 
		}
		std::vector<int> MC_move{a,b,c};
	   	vec_neigh[MC_move]=d;
	}

	return vec_neigh;
}


std::map<std::pair<int,int>,double> neigh_angles=initialize_map_neigh_angles();

std::map<std::vector<int>,int> MC_move_n = initialize_map_vec_neigh();


#endif
