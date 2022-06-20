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


		polymer(int &l,int &n):length{l},number_of_polymer{n}{}
		
		polymer(int &l,int &n,int &t,int &u,int &v):length{l},number_of_polymer{n}{ 
			for(auto i=0;i<length;i++){
				auto tmp=monomer(t,u,v);
				this->chain.push_back(tmp);
			}
		}
		
		polymer(const polymer& p){
			length=p.length;
			for(auto i=0;i<length;i++){
				chain.push_back(p.chain[i]);
			}
		}

		void operator=(polymer& p){
			for(auto i=0;i<p.chain.size();i++){
				this->chain.push_back(p.chain[i]);
			}
		}	
}; 


#endif