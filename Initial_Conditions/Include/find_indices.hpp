#ifndef __find_indices__
#define __find_indices__
#include "observables.hpp"

int find_neighrest_neighbors_left(polymer& p,int index_monomer){
	int nn_left_index=index_monomer-1;
	while(nn_left_index>=0){
		double d= distance(p.chain[index_monomer],p.chain[nn_left_index]);
		if(fabs(d-1)<0.001){
			return nn_left_index;
		}
		else{
			nn_left_index-=1;
		}
	}
	return -1;
}

int find_neighrest_neighbors_right(polymer& p,int index_monomer){
	int nn_right_index=index_monomer+1;
	if(index_monomer==-1){
		return -1;
	}
	else{
		while(nn_right_index<=int(p.chain.size())-1){
			double d = distance(p.chain[index_monomer],p.chain[nn_right_index]);
			if(fabs(d-1)<0.001){
				return nn_right_index;
			}
			else{
				nn_right_index+=1;
			}
		}
		return -1;
	}
}

#endif