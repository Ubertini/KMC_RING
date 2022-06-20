#include <string>
#include <vector>
#include "starting_points.hpp"
#include "monte_carlo_real_ring.hpp"

double bending_rigidity = 2.;

int main(){
	
	std::string type_polymer = "Real";

	int polymer_length=40;
	int number_of_polymers=1000;
	double monomer_density=1.25;
	std::string density= std::to_string(int(monomer_density*100));
	
	long long int MC_steps=1e6;           //6e9; //2e9;  //1e9;
	int stride = 1e5;
	int number_of_runs=1;
	int start_run=0;
	
	bool first_run=1;

	std::string DIR_MASTER = "DATA";
	std::string command0 = "mkdir -p "+DIR_MASTER;
	system(command0.c_str());

	std::vector<int> start_traj_vec{};

	for(int run=start_run;run<number_of_runs;run++){
		std::string DIR = DIR_MASTER+"/"+std::to_string(run)+"/";
		std::cout<<DIR<<std::endl;
		std::string command = "mkdir -p "+DIR;
		system(command.c_str());
		
		if(first_run){
			start_traj_vec.push_back(0);	
		}
		else{
			int N_possible_step=300000;
			start_traj_vec.push_back(find_start_points(DIR,type_polymer,polymer_length,density,
			MC_steps,stride,run,N_possible_step));
		}
	}
	
	if(start_traj_vec.size()!=number_of_runs){
		return 1;
	}

	else{
		MC_routine_real_ring(DIR_MASTER,polymer_length,number_of_polymers,monomer_density,start_traj_vec,MC_steps,stride,start_run,number_of_runs);
	}
	
	return 0;
}
