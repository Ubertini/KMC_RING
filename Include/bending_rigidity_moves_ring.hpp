#ifndef __bending_rigidity_moves_ring__
#define __bending_rigidity_moves_ring__

#include "observables.hpp"
#include "initialize_maps.hpp"
#include "initialize_MC_moves.hpp"
#include <math.h>   

extern std::map<std::pair<int,int>,double> neigh_angles;
extern std::map<std::vector<int>,int> MC_move_n;
extern double bending_rigidity;

bool Rouse_move_central_monomers_from_left_stored_length_ring(monomer& trial,polymer &p,const int &index_monomer,const int &index_right,const int &index_left,double &E){
	
	int r=index_right;
	int rr;
	int l;

	//FIND rr
	if(r<int(p.chain.size())-1){
		if(p.chain[r+1].stored_length_left){
			if(r<int(p.chain.size())-2)
				rr=r+2;
			else
				rr=0;
		}
		else
			rr=r+1;	
	}
	else{
		if(p.chain[0].stored_length_left)
			rr=1;
		else
			rr=0;
	}

	//FIND l
	if(index_left==0)
		l=int(p.chain.size())-1;
	else
		l=index_left-1;	

	std::vector<int> t_r = p.chain[index_monomer]-p.chain[r];
	int MC_tr=MC_move_n[t_r];

	double MH=0;
	double deltaE=0;

	std::vector<int> t_rr = p.chain[r]-p.chain[rr];
	int MC_trr=MC_move_n[t_rr];
	std::vector<int> t_l = p.chain[l]-p.chain[index_monomer];
	int MC_tl=MC_move_n[t_l];
	auto pair_2=std::make_pair(MC_tl,MC_tr);			
	auto pair_3=std::make_pair(MC_tr,MC_trr);
	double cos_theta_initial= neigh_angles[pair_2]+neigh_angles[pair_3];
	
	std::vector<int> t_l_new = p.chain[index_monomer]-trial;					
	int MC_tl_new=MC_move_n[t_l_new];
	std::vector<int> t_r_new = trial-p.chain[r];					
	int MC_tr_new=MC_move_n[t_r_new];
	auto pair_1_new=std::make_pair(MC_tl,MC_tl_new);
	auto pair_2_new=std::make_pair(MC_tl_new,MC_tr_new);
	auto pair_3_new=std::make_pair(MC_tr_new,MC_trr);
	double cos_theta_final = neigh_angles[pair_1_new]+neigh_angles[pair_2_new]+neigh_angles[pair_3_new];
	
	deltaE = bending_rigidity*(1+cos_theta_initial-cos_theta_final);
	MH = exp(-deltaE);		
	
	double random_numer = ((double) rand() / (RAND_MAX));

	if(random_numer<=MH){
		E+=deltaE;
		p.chain[index_monomer].stored_length_left=0;
		p.chain[index_left].stored_length_right=0;
		return true;
	}
	else
		return false;
}

bool Rouse_move_central_monomers_from_right_stored_length_ring(monomer& trial,polymer &p,const int &index_monomer,const int &index_right,const int &index_left,double &E){
	
	int r;
	int l=index_left;
	int ll;

	//FIND r
	if(index_right==int(p.chain.size())-1)
		r=0;
	else
		r=index_right+1;	
	
	//FIND ll
	if(l>0){
		if(p.chain[l-1].stored_length_right){
			if(l>1)
				ll=l-2;
			else
				ll=int(p.chain.size())-1;
		}
		else
			ll=l-1;	
	}
	else{
		if(p.chain[int(p.chain.size())-1].stored_length_right)
			ll=int(p.chain.size())-2;
		else
			ll=int(p.chain.size())-1;
	}
							
	std::vector<int> t_l = p.chain[l]-p.chain[index_monomer];
	int MC_tl=MC_move_n[t_l];

	double MH=0;
	double deltaE=0;

	std::vector<int> t_ll = p.chain[ll]-p.chain[l];
	int MC_tll=MC_move_n[t_ll];
	std::vector<int> t_r = p.chain[index_monomer]-p.chain[r];
	int MC_tr=MC_move_n[t_r];
	auto pair_1=std::make_pair(MC_tll,MC_tl);			
	auto pair_2=std::make_pair(MC_tl,MC_tr);
	double cos_theta_initial= neigh_angles[pair_1]+neigh_angles[pair_2];
	
	std::vector<int> t_l_new = p.chain[l]-trial;					
	int MC_tl_new=MC_move_n[t_l_new];
	std::vector<int> t_r_new = trial-p.chain[index_monomer];					
	int MC_tr_new=MC_move_n[t_r_new];
	auto pair_1_new=std::make_pair(MC_tll,MC_tl_new);
	auto pair_2_new=std::make_pair(MC_tl_new,MC_tr_new);
	auto pair_3_new=std::make_pair(MC_tr_new,MC_tr);
	double cos_theta_final = neigh_angles[pair_1_new]+neigh_angles[pair_2_new]+neigh_angles[pair_3_new];
	
	deltaE = bending_rigidity*(1+cos_theta_initial-cos_theta_final);
	MH = exp(-deltaE);

	double random_numer = ((double) rand() / (RAND_MAX));

	if(random_numer<=MH){
		E+=deltaE;
		p.chain[index_monomer].stored_length_right=0;
		p.chain[index_right].stored_length_left=0;
		return true;
	}
	else
		return false;
}

//UPDATED WITH STORED LENGTHS
bool Rouse_move_central_monomers_no_stored_length_ring(monomer& trial,polymer &p,const int &index_monomer,const int &index_right,const int &index_left,double &E){
	
	int r=index_right;
	int rr;
	int l=index_left;
	int ll;

	//FIND rr
	if(r<int(p.chain.size())-1){
		if(p.chain[r+1].stored_length_left){
			if(r<int(p.chain.size())-2)
				rr=r+2;
			else
				rr=0;
		}
		else
			rr=r+1;	
	}
	else{
		if(p.chain[0].stored_length_left)
			rr=1;
		else
			rr=0;
	}
	
	//FIND ll
	if(l>0){
		if(p.chain[l-1].stored_length_right){
			if(l>1)
				ll=l-2;
			else
				ll=int(p.chain.size())-1;
		}
		else
			ll=l-1;	
	}
	else{
		if(p.chain[int(p.chain.size())-1].stored_length_right)
			ll=int(p.chain.size())-2;
		else
			ll=int(p.chain.size())-1;
	}

	std::vector<int> t_l = p.chain[l]-p.chain[index_monomer];
	int MC_tl=MC_move_n[t_l];
	std::vector<int> t_r = p.chain[index_monomer]-p.chain[r];
	int MC_tr=MC_move_n[t_r];
	
	double MH=0;
	double deltaE=0;

	std::vector<int> t_ll = p.chain[ll]-p.chain[l];
	int MC_tll=MC_move_n[t_ll];
	std::vector<int> t_rr = p.chain[r]-p.chain[rr];
	int MC_trr=MC_move_n[t_rr];
	auto pair_1=std::make_pair(MC_tll,MC_tl);
	auto pair_2=std::make_pair(MC_tl,MC_tr);
	auto pair_3=std::make_pair(MC_tr,MC_trr);
	double cos_theta_initial= neigh_angles[pair_1]+neigh_angles[pair_2]+neigh_angles[pair_3];
	
	std::vector<int> t_l_new = p.chain[l]-trial;					
	int MC_tl_new=MC_move_n[t_l_new];
	std::vector<int> t_r_new = trial-p.chain[r];					
	int MC_tr_new=MC_move_n[t_r_new];
	auto pair_1_new=std::make_pair(MC_tll,MC_tl_new);
	auto pair_2_new=std::make_pair(MC_tl_new,MC_tr_new);
	auto pair_3_new=std::make_pair(MC_tr_new,MC_trr);
	double cos_theta_final = neigh_angles[pair_1_new]+neigh_angles[pair_2_new]+neigh_angles[pair_3_new];
	
	deltaE = bending_rigidity*(cos_theta_initial-cos_theta_final);
	MH = exp(-deltaE);
	
	double random_numer = ((double) rand() / (RAND_MAX));

	if(random_numer<=MH){
		E+=deltaE;
		return true;
	}
	else
		return false;
}

bool retraction_move_central_monomers_ring(polymer &p,const int &index_monomer,const int &index_right,const int &index_left,double &E){
	
	int r=index_right;
	int rr;
	int l=index_left;
	int ll;

	//FIND rr
	if(r<int(p.chain.size())-1){
		if(p.chain[r+1].stored_length_left){
			if(r<int(p.chain.size())-2)
				rr=r+2;
			else
				rr=0;
		}
		else
			rr=r+1;	
	}
	else{
		if(p.chain[0].stored_length_left)
			rr=1;
		else
			rr=0;
	}
	
	//FIND ll
	if(l>0){
		if(p.chain[l-1].stored_length_right){
			if(l>1)
				ll=l-2;
			else
				ll=int(p.chain.size())-1;
		}
		else
			ll=l-1;	
	}
	else{
		if(p.chain[int(p.chain.size())-1].stored_length_right)
			ll=int(p.chain.size())-2;
		else
			ll=int(p.chain.size())-1;
	}
							
	std::vector<int> t_l = p.chain[l]-p.chain[index_monomer];
	int MC_tl=MC_move_n[t_l];
	std::vector<int> t_r = p.chain[index_monomer]-p.chain[r];
	int MC_tr=MC_move_n[t_r];
	
	double MH=0;
	double deltaE=0;

	std::vector<int> t_ll = p.chain[ll]-p.chain[l];
	int MC_tll=MC_move_n[t_ll];
	std::vector<int> t_rr = p.chain[r]-p.chain[rr];
	int MC_trr=MC_move_n[t_rr];
	auto pair_1=std::make_pair(MC_tll,MC_tl);
	auto pair_2=std::make_pair(MC_tl,MC_tr);
	auto pair_3=std::make_pair(MC_tr,MC_trr);
	double cos_theta_initial= neigh_angles[pair_1]+neigh_angles[pair_2]+neigh_angles[pair_3];
	
	std::vector<int> t_new = p.chain[l]-p.chain[r];					
	int MC_new=MC_move_n[t_new];
	auto pair_1_new=std::make_pair(MC_tll,MC_new);
	auto pair_2_new=std::make_pair(MC_new,MC_trr);
	double cos_theta_final = neigh_angles[pair_1_new]+neigh_angles[pair_2_new];
	deltaE = bending_rigidity*(-1+cos_theta_initial-cos_theta_final);
	MH = exp(-deltaE);

	double random_numer = ((double) rand() / (RAND_MAX));

	if(random_numer<=MH){
		E+=deltaE;
		return true;
	}
	else
		return false;
}
#endif
