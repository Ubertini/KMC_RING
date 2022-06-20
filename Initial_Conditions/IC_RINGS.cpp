#include <fstream>
#include <string>
#include <iostream>
#include "polymer.hpp"
#include "initial_condition.hpp"

int main(){
	std::cout<<"Lattice size = "<<Lattice_size<<std::endl;
	std::string type_polymer = "Ring";
	std::vector<int> number_of_chain{40};
	std::vector<int> N{1000};

	for(int i=0;i<N.size();i+=1){
		std::cout<<"Number of monomers = "<<N[i]<<" number of chain = "<<number_of_chain[i]<<std::endl;
		if(type_polymer=="Ring"){
			Initial_Condition_Ring(N[i],number_of_chain[i]);
		}
		else{
			Initial_Condition(N[i],number_of_chain[i]);
		}
	}	
	return 0;
}
