#ifndef __MC_moves__
#define __MC_moves__
#include <string>
#include <sstream>

struct MC_move{		
	int dt;
	int du;
	int dv;
	MC_move(int t,int u,int v):dt{t},du{u},dv{v} {};
};

std::vector<MC_move> initialize_MC_move(){
	
	std::vector<MC_move> mc_move;
	std::ifstream infile("Infile/MC_moves.txt");		

	std::string line;
	while (std::getline(infile, line))
	{
	    std::istringstream iss(line);
	    int a, b, c, d;
	    if (!(iss >> a >> b >> c>>d)) { 
	    	break; 
		}	 	
	mc_move.push_back(MC_move(a,b,c));  	
	}
	return mc_move;
}

#endif