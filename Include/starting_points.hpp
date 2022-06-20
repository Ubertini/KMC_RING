#ifndef __starting_points__
#define __starting_points__

#include <fstream>
#include <string>
#include <iostream>

extern double bending_rigidity;

bool is_file_exist(std::string fileName)
{
    std::ifstream infile(fileName.c_str());
    return infile.good();
}

std::vector<int> find_start_points(std::string type_polymer,int polymer_length,std::string density,
	long long int MC_steps,int stride,int start_run,int number_of_runs,int N_possible_step){
	
	std::string DIR_MASTER = "/net/sbp/sbpstore1/mubertin/DATA_CUBIC/RINGS/Real_"+density+"_"
	+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length);
	
	std::vector<int> start_traj_vec;

	for(int i=0;i<number_of_runs;i++){
		std::string DIR = DIR_MASTER+"/"+std::to_string(i+start_run)+"/";
		int start_traj_temp=-1;
		for(int step=0;step<N_possible_step;step++){
			std::string fileName = DIR+"traj_"+"Real"+"_"+density+"_"
			+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)
			+"_"+std::to_string(step)+"_"+std::to_string(i+start_run);
			if(is_file_exist(fileName)){
				start_traj_temp=step;	
			}
			else{
				break;
			}
		}
		if(start_traj_temp == -1){
			std::cout<<"I never simulated a system at this condition hence I do not find any configuration"<<std::endl;
			return start_traj_vec;
		}
		else{
			start_traj_vec.push_back(start_traj_temp);
		}
	}
	return start_traj_vec;
}

#endif