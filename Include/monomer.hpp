#ifndef __monomer__
#define __monomer__
#include <cmath>
#include <vector>
#include "lattice.hpp"


int minimum_image_coordinate(int coordinate){
	if(coordinate >= Lattice_size){
		while(coordinate >= Lattice_size){
			coordinate-= Lattice_size;
		}
	}
	else{
		while(coordinate < 0){
			coordinate += Lattice_size;
		}
	} 
	return coordinate;
}

struct monomer{
	int t; //three unit vector of the fcc lattice which form our basis
	int u;
	int v;

	bool stored_length_left;
	bool stored_length_right;

	int t_image;
	int u_image;
	int v_image;	

	monomer(const int& a,const int& b,const int& c): t{a},u{b},v{c}{
		stored_length_right=0;
		stored_length_left=0;
		if(t >= Lattice_size || t < 0)
			t_image = minimum_image_coordinate(t);
		else
			t_image = t;
		
		if(u >= Lattice_size || u < 0)
			u_image=minimum_image_coordinate(u);
		else
			u_image = u;
		if(v >= Lattice_size || v < 0)
			v_image=minimum_image_coordinate(v);
		else
			v_image = v;
	};

	void operator=(monomer& m){
		t=m.t;
		u=m.u;
		v=m.v;
		stored_length_left=m.stored_length_left;
		stored_length_right=m.stored_length_right;
		t_image=m.t_image;
		u_image=m.u_image;
		v_image=m.v_image;

		
	};

	std::vector<int> operator-(monomer& m){
		int x =t-m.t;
		int y = u-m.u;
		int z = v-m.v;
		std::vector<int> diff{x,y,z};
		return diff;
	};

	double cartesian_coordinate_x(){
		return 0.70710678118*(t);
	}
	double cartesian_coordinate_y(){
		return 0.70710678118*(u);
	}
	double cartesian_coordinate_z(){
		return 0.70710678118*(v);
	}

	double cartesian_coordinate_x_image(){
		return 0.70710678118*(t_image);
	}
	double cartesian_coordinate_y_image(){
		return 0.70710678118*(u_image);
	}
	double cartesian_coordinate_z_image(){
		return 0.70710678118*(v_image);
	}		

};


double distance(monomer& a, monomer& b){
		double distance= ((a.cartesian_coordinate_x()-b.cartesian_coordinate_x())*(a.cartesian_coordinate_x()-b.cartesian_coordinate_x())
			+(a.cartesian_coordinate_y()-b.cartesian_coordinate_y())*(a.cartesian_coordinate_y()-b.cartesian_coordinate_y())
			+(a.cartesian_coordinate_z()-b.cartesian_coordinate_z())*(a.cartesian_coordinate_z()-b.cartesian_coordinate_z()));
		return distance;
}

double distance_images(monomer& a, monomer& b){
		double distance=((a.cartesian_coordinate_x_image()-b.cartesian_coordinate_x_image())*(a.cartesian_coordinate_x_image()-b.cartesian_coordinate_x_image())
			+(a.cartesian_coordinate_y_image()-b.cartesian_coordinate_y_image())*(a.cartesian_coordinate_y_image()-b.cartesian_coordinate_y_image())
			+(a.cartesian_coordinate_z_image()-b.cartesian_coordinate_z_image())*(a.cartesian_coordinate_z_image()-b.cartesian_coordinate_z_image()));
		return distance;
}










#endif
