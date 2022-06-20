#include <string>
#include <vector>
#include "starting_points.hpp"
#include "monte_carlo_real_ring.hpp"

double bending_rigidity = 1.;

int main(){
	
	std::string type_polymer = "Real";

	int polymer_length=1280;
	int number_of_polymers=31;
	double monomer_density=1.25;
	std::string density= std::to_string(int(monomer_density*100));
	
	long long int MC_steps=3e13;           //6e9; //2e9;  //1e9;
	int stride = 1e9;
	int number_of_runs=1;
	int start_run=0;
	
	bool first_run=0;
	
	std::vector<int> start_traj_vec{};
	if(first_run){
		for(int i=0;i<number_of_runs;i++){
			start_traj_vec.push_back(0);
		}
	}
	else{
		int N_possible_step=300000;
		start_traj_vec=find_start_points(type_polymer,polymer_length,density,
		MC_steps,stride,start_run,number_of_runs,N_possible_step);
	}
	
	if(start_traj_vec.size()!=number_of_runs){
		return 1;
	}
	else{
		MC_routine_real_ring(polymer_length,number_of_polymers,monomer_density,start_traj_vec,MC_steps,stride,start_run,number_of_runs);
	}
	
	return 0;
}
