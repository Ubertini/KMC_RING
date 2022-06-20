#ifndef __observables__
#define __observables__
#include "polymer.hpp"
#include "initialize_MC_moves.hpp"
#include "initialize_maps.hpp"

extern double bending_rigidity;

/*
double compute_bending_energy(polymer& p){
	int r_old=-1;
	int l_old=-1;
	double E =0;
	for(int j=1;j<p.chain.size();j++){
		int l=find_neighrest_neighbors_left(p,j);
		int r=find_neighrest_neighbors_right(p,j);
		if(l!=-1 && r !=-1){
			if(r!=r_old && l!=l_old){
				std::vector<int> t_l = p.chain[l]-p.chain[j];
				int MC_tl=MC_move_n[t_l];
				std::vector<int> t_r = p.chain[j]-p.chain[r];
				int MC_tr=MC_move_n[t_r];
				auto pair_2=std::make_pair(MC_tl,MC_tr);
				double cos_theta=neigh_angles[pair_2];
				E+=bending_rigidity*(1-cos_theta);
				r_old=r;
				l_old=l;
			}
		}
	}
	return E;
}
*/

double gyration_radius(polymer a){
	//compute the centre of mass
	std::vector<double> v{0,0,0};
	for(unsigned int  i=0;i<a.chain.size();i++){
		v[0]+=a.chain[i].cartesian_coordinate_x();
		v[1]+=a.chain[i].cartesian_coordinate_y();
		v[2]+=a.chain[i].cartesian_coordinate_z();
	}
	v[0]/=a.chain.size();
	v[1]/=a.chain.size();
	v[2]/=a.chain.size();
	
	double gyration_radius=0;

	for(unsigned int  i=0;i<a.chain.size();i++){
		gyration_radius+=(a.chain[i].cartesian_coordinate_x()-v[0])*(a.chain[i].cartesian_coordinate_x()-v[0])
						+(a.chain[i].cartesian_coordinate_y()-v[1])*(a.chain[i].cartesian_coordinate_y()-v[1])
						+(a.chain[i].cartesian_coordinate_z()-v[2])*(a.chain[i].cartesian_coordinate_z()-v[2]);
	}

	gyration_radius/=double(a.chain.size());
	return gyration_radius;
}

double end_end_distance(polymer a){
	return distance(a.chain[0],a.chain.back());
}

double G1(polymer reference, polymer p){
	double g1;
	
	for(std::size_t i=0;i<reference.chain.size();i++){
		g1+=distance(reference.chain[i],p.chain[i]);
	}

	return g1/=double(reference.chain.size());
}

double G2(polymer reference, polymer p,std::vector<double> r_cm_0){
	double g2;

	std::vector<double> r_cm_t{0,0,0};
	for(unsigned int  i=0;i<p.chain.size();i++){
		r_cm_t[0]+=p.chain[i].cartesian_coordinate_x();
		r_cm_t[1]+=p.chain[i].cartesian_coordinate_y();
		r_cm_t[2]+=p.chain[i].cartesian_coordinate_z();
	}
	
	r_cm_t[0]/=p.chain.size();
	r_cm_t[1]/=p.chain.size();
	r_cm_t[2]/=p.chain.size();

	for(std::size_t i=0;i<reference.chain.size();i++){
		g2+= (p.chain[i].cartesian_coordinate_x()-r_cm_t[0]-reference.chain[i].cartesian_coordinate_x()+r_cm_0[0])*(p.chain[i].cartesian_coordinate_x()-r_cm_t[0]-reference.chain[i].cartesian_coordinate_x()+r_cm_0[0])
			+ (p.chain[i].cartesian_coordinate_y()-r_cm_t[1]-reference.chain[i].cartesian_coordinate_y()+r_cm_0[1])*(p.chain[i].cartesian_coordinate_y()-r_cm_t[1]-reference.chain[i].cartesian_coordinate_y()+r_cm_0[1])+
			+(p.chain[i].cartesian_coordinate_z()-r_cm_t[2]-reference.chain[i].cartesian_coordinate_z()+r_cm_0[2])*(p.chain[i].cartesian_coordinate_z()-r_cm_t[2]-reference.chain[i].cartesian_coordinate_z()+r_cm_0[2]);
	}
	return g2/=double(reference.chain.size());
}

double G3(std::vector<double> r_cm_0, polymer p){
	double g3;

	std::vector<double> r_cm_t{0,0,0};
	for(unsigned int  i=0;i<p.chain.size();i++){
		r_cm_t[0]+=p.chain[i].cartesian_coordinate_x();
		r_cm_t[1]+=p.chain[i].cartesian_coordinate_y();
		r_cm_t[2]+=p.chain[i].cartesian_coordinate_z();
	}
	
	r_cm_t[0]/=p.chain.size();
	r_cm_t[1]/=p.chain.size();
	r_cm_t[2]/=p.chain.size();

	for(std::size_t i=0;i<p.chain.size();i++){
		g3+= (r_cm_t[0]-r_cm_0[0])*(r_cm_t[0]-r_cm_0[0])+(r_cm_t[1]-r_cm_0[1])*(r_cm_t[1]-r_cm_0[1])+
			+(r_cm_t[2]-r_cm_0[2])*(r_cm_t[2]-r_cm_0[2]);
	}
	return g3/=double(p.chain.size());
}

std::vector<double> compute_R_cm(polymer p){
	std::vector<double> r_cm_t{0,0,0};
	for(unsigned int  i=0;i<p.chain.size();i++){
		r_cm_t[0]+=p.chain[i].cartesian_coordinate_x();
		r_cm_t[1]+=p.chain[i].cartesian_coordinate_y();
		r_cm_t[2]+=p.chain[i].cartesian_coordinate_z();
	}
	r_cm_t[0]/=p.chain.size();
	r_cm_t[1]/=p.chain.size();
	r_cm_t[2]/=p.chain.size();
	return r_cm_t;
}

#endif
