#ifndef __acceptance_ring__
#define __acceptance_ring__
#include "observables.hpp"
#include "initialize_maps.hpp"
#include "initialize_MC_moves.hpp"
#include "bending_rigidity_moves_ring.hpp"
#include <math.h>   

extern std::map<std::pair<int,int>,double> neigh_angles;
extern std::map<std::vector<int>,int> MC_move_n;
extern double bending_rigidity;

inline bool distance_real_ring(int (&lattice)[Lattice_size][Lattice_size][Lattice_size],monomer& trial, polymer& p, unsigned int index_monomer,double& E){
		
		double distance_plus;
		double distance_minus;
		int index_left;
		int index_right;

		if(index_monomer>0 && index_monomer < int(p.chain.size())-1){
			index_left=index_monomer-1;
			index_right=index_monomer+1;
			distance_plus = distance(trial,p.chain[index_right]);
			distance_minus = distance(trial,p.chain[index_left]);
		}

		else if(index_monomer==0){
			index_left=int(p.chain.size())-1;
			index_right=1;
			//std::cout<<index_left<<"\t"<<index_right<<std::endl;
			distance_plus = distance(trial,p.chain[1]);
			distance_minus = distance(trial,p.chain[index_left]);
			//std::cout<<"Plus = "<<distance_plus<<std::endl;
			//std::cout<<"Minus = "<<distance_minus<<std::endl;
		}

		else{
			index_left=index_monomer-1;
			index_right=0;
			distance_plus =  distance(trial,p.chain[0]);
			distance_minus = distance(trial,p.chain[index_left]);
			//std::cout<<distance_plus<<std::endl;
		}

		if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus-1)<0.0001 
			&& lattice[trial.t_image][trial.u_image][trial.v_image]==0){
			//Rolling out a left stored length
			if(p.chain[index_monomer].stored_length_left and !(p.chain[index_monomer].stored_length_right)){
				return Rouse_move_central_monomers_from_left_stored_length_ring(trial,p,index_monomer,index_right,index_left,E);
			}
			else if((!p.chain[index_monomer].stored_length_left) and p.chain[index_monomer].stored_length_right){
				return Rouse_move_central_monomers_from_right_stored_length_ring(trial,p,index_monomer,index_right,index_left,E);
			}
			else{ //Rouse move from a non stored length
				return Rouse_move_central_monomers_no_stored_length_ring(trial,p,index_monomer,index_right,index_left,E);
			}
		}

		else if(fabs(distance_minus)<0.0001 && fabs(distance_plus-1)<0.0001){
			if(!(p.chain[index_left].stored_length_left)){
				if(p.chain[index_monomer].stored_length_right){ //diffusion along the backbone we accept it
					p.chain[index_monomer].stored_length_right=0; //towards left
					p.chain[index_right].stored_length_left=0;
					p.chain[index_monomer].stored_length_left=1;
					p.chain[index_left].stored_length_right=1;
					return true; 
				}
				else{ // -> retraction move we have to check angles 
					bool acceptance=retraction_move_central_monomers_ring(p,index_monomer,index_right,index_left,E);
					if(acceptance){
						p.chain[index_monomer].stored_length_left=1;
						p.chain[index_left].stored_length_right=1;
					}
					return acceptance;
				}
			}
			else
				return false;
		}
		
		else if(fabs(distance_plus)<0.0001 && fabs(distance_minus-1)<0.0001){
			
			if(!(p.chain[index_right].stored_length_right)){
					if(p.chain[index_monomer].stored_length_left){ //diffusion along the backbone
						p.chain[index_monomer].stored_length_left=0;
						p.chain[index_left].stored_length_right=0;
						p.chain[index_monomer].stored_length_right=1;
						p.chain[index_right].stored_length_left=1;
						return true; 
					}
					else{ // -> retraction move we have to check angles 
						bool acceptance=retraction_move_central_monomers_ring(p,index_monomer,index_right,index_left,E);
						if(acceptance){
							p.chain[index_monomer].stored_length_right=1;
							p.chain[index_right].stored_length_left=1;
						}
						return acceptance;
					}
			}

			else 
				return false;
		}
		
		else
			return false;
	
}

inline bool acceptance_real_ring(int (&lattice)[Lattice_size][Lattice_size][Lattice_size],
	monomer trial,polymer& p,const unsigned int index_monomer,double& E){
		return distance_real_ring(lattice,trial,p,index_monomer,E);
}



#endif
