#ifndef __polymer__
#define __polymer__

#include "monomer.hpp"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <vector>
#include <fstream>
#include <iostream>

struct polymer{
	
	int length;
	int number_of_polymer;
	std::vector<monomer> chain;
	std::vector<double> angles;

	//constructor setting all monomers in the same random lattice site
	/*polymer(int &l):length{l}{ 
		srand(time(NULL));
		int t= rand() % Lattice_size;
		int u= rand() % Lattice_size;
		int v= rand() % Lattice_size;
		for(auto i=0;i<length;i++){
			auto tmp=monomer(t,u,v);
			this->chain.push_back(tmp);
			}
	}*/
		polymer(int &l,int &n):length{l},number_of_polymer{n}{}
		
		polymer(int &l,int &n,int &t,int &u,int &v):length{l},number_of_polymer{n}{ 
			
			for(auto i=0;i<length;i++){
				
				auto tmp=monomer(t,u,v);
				this->chain.push_back(tmp);
			
			}
			
			for(auto i=0;i<length;i++){
				
			}
		}


		polymer(const polymer& p){
			length=p.length;
			for(auto i=0;i<length;i++){
				chain.push_back(monomer(p.chain[i].t,p.chain[i].u,p.chain[i].v));
			}
		}	
}; 


#endif