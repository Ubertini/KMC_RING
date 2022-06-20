#ifndef __relaxation__
#define __relaxation__
#include <string>
#include "acceptance_linear.hpp"
#include "initialize_MC_moves.hpp"
#include "initialize_maps.hpp"
#include <random>

void relaxation(std::vector<polymer> &melt,int (&lattice)[Lattice_size][Lattice_size][Lattice_size],long long int MC_steps,double& E){
	std::vector<MC_move> mc_move=initialize_MC_move();
	for(long long int j=0;j<MC_steps;j++){ 				
			
			unsigned int polymer_to_move = rand()%melt.size();
			unsigned int monomer_to_move = rand()%melt[polymer_to_move].chain.size();
			auto MC_move=rand()%12; //there are 12 possible movement
			monomer trial{melt[polymer_to_move].chain[monomer_to_move].t+mc_move[MC_move].dt,melt[polymer_to_move].chain[monomer_to_move].u+mc_move[MC_move].du,melt[polymer_to_move].chain[monomer_to_move].v+mc_move[MC_move].dv};
				
			if(acceptance_ideal_linear(trial,melt[polymer_to_move],monomer_to_move,E)){
					
				lattice[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]-=1;
				melt[polymer_to_move].chain[monomer_to_move].t=trial.t;
				melt[polymer_to_move].chain[monomer_to_move].u=trial.u;
				melt[polymer_to_move].chain[monomer_to_move].v=trial.v;
				melt[polymer_to_move].chain[monomer_to_move].t_image=trial.t_image;
				melt[polymer_to_move].chain[monomer_to_move].u_image=trial.u_image;
				melt[polymer_to_move].chain[monomer_to_move].v_image=trial.v_image;
				lattice[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]+=1;
					
				/*double E_check=0;
				for(int n_pol=0;n_pol<melt.size();n_pol++){
					E_check+=compute_bending_energy(melt[n_pol]);
				}
					
				if(E_check!=E){
					std::cout<<"ERROR energy computed now vs Energy updated: "<<E_check<<"\t"<<E<<"\n";
				}*/
		}			
	}
}

#endif