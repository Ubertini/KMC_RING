#ifndef __monte_carlo_real_ring__
#define __monte_carlo_real_ring__
#include <string>
#include <random>
#include <dirent.h>
#include <errno.h>
#include <fstream>
#include <string>
#include <iostream>
#include "timer.hpp"
timer<> t;
#include "acceptance_ring.hpp"
#include "initialize_MC_moves.hpp"
#include "initialize_maps.hpp"

extern std::map<std::pair<int,int>,double> neigh_angles;
extern std::map<std::vector<int>,int> MC_move_n;
extern double bending_rigidity;

void MC_routine_real_ring(std::string DIR_MASTER,int polymer_length, int number_of_polymers,
	double monomer_density,std::vector<int> start_traj_vec, long long int MC_steps, int stride,int start_run,int number_of_runs){
	
	std::string density;

	if(number_of_polymers==1){
		density="0";
	}
	
	else{
		density= std::to_string(int(monomer_density*100));
	}
	
	std::vector<MC_move> mc_move=initialize_MC_move();

	for(auto run=start_run;run<start_run+number_of_runs;run++){

		std::string DIR = DIR_MASTER+"/"+std::to_string(run)+"/";
		
		int start_traj=start_traj_vec[run-start_run];
		
		std::vector<polymer> melt;
		
		int lattice[Lattice_size][Lattice_size][Lattice_size]={0};
		double E = 0; // Energy associated to bending rigidity
		
		std::ifstream file;

		if(start_traj==0){
			if(number_of_polymers==1){
				std::cout<<"ICS_RINGS/IC_RING_"+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)<<std::endl; 
				file.open("ICS_RINGS/IC_RING_"+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)); 
			}
			else{
				std::cout<<"ICS_RINGS/IC_RING_"+density+"_"+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)<<std::endl; 
				file.open("ICS_RINGS/IC_RING_"+density+"_"+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)); 
			}	
		}
		
		else{
			std::cout<<DIR+"traj_Real_"+density+"_"+std::to_string(int(bending_rigidity))
			+"_"+std::to_string(polymer_length)+"_"
			+std::to_string(start_traj)+"_"+std::to_string(run)<<std::endl;
			file.open(DIR+"traj_Real_"+density+"_"+std::to_string(int(bending_rigidity))
			+"_"+std::to_string(polymer_length)+"_"
			+std::to_string(start_traj)+"_"+std::to_string(run));
		}
		
		std::vector<int> coordinates;
		double val=0;
		if(file.is_open()){
			while(!file.eof()){
	      		while(file >> val){
	        		coordinates.push_back(int(val));
	      		}
	    	}
	  	}
	  	else{
	       std::cout << "unable to open file."<<std::endl;
	  	}
	  	file.close();

	  	for(int n_pol=0;n_pol<number_of_polymers;n_pol++){
			polymer polymer_tmp{polymer_length,n_pol};
			for(int j=0;j<polymer_length;j++){
				int t = coordinates[3*polymer_length*n_pol+3*j];
				int u = coordinates[3*polymer_length*n_pol+3*j+1];
				int v = coordinates[3*polymer_length*n_pol+3*j+2];
				monomer monomer_tmp{t,u,v};
				polymer_tmp.chain.push_back(monomer_tmp);
				if(j>0){
					double d = distance(polymer_tmp.chain[j],polymer_tmp.chain[j-1]);
					auto t_l = polymer_tmp.chain[j]-polymer_tmp.chain[j-1];
					int diff_t = polymer_tmp.chain[j].t-polymer_tmp.chain[j-1].t-t_l[0];
					int diff_u = polymer_tmp.chain[j].u-polymer_tmp.chain[j-1].u-t_l[1];
					int diff_v = polymer_tmp.chain[j].v-polymer_tmp.chain[j-1].v-t_l[2];
					
					if(diff_t!=0 || diff_u!=0 || diff_v!=0){
						std::cout<<diff_t<<"\t"<<diff_u<<"\t"<<diff_v<<std::endl;
					}

					if(fabs(d)<0.001){
						polymer_tmp.chain[j].stored_length_left=1;
						polymer_tmp.chain[j-1].stored_length_right=1;
					}
				}
				lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]+=1;
			}
			double d = distance(polymer_tmp.chain[0],polymer_tmp.chain[polymer_length-1]);
			if(fabs(d)<0.001){
				polymer_tmp.chain[0].stored_length_left=1;
				polymer_tmp.chain[polymer_length-1].stored_length_right=1;
			}
			
			melt.push_back(polymer_tmp);
		}

		srand(time(NULL));

		std::string name_file0= DIR+"E_bending_"+"Real"+"_"+density+"_"
		+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)+"_"+std::to_string(run);
		std::ofstream E_bending(name_file0);


		std::string name_file1= DIR+"GR_"+"Real"+"_"+density+"_"
		+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)+"_"+std::to_string(run);
		std::ofstream GR(name_file1);

		std::string name_file2= DIR+"R_EE_"+"Real"+"_"+density+"_"
		+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)+"_"+std::to_string(run);
		std::ofstream R_EE(name_file2);

		std::string name_file3 = DIR+"traj_"+"Real"+"_"+density+"_"
		+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)+"_"
		+std::to_string(start_traj)+"_"+std::to_string(run);
		std::ofstream traj;
		traj.open(name_file3);

		std::string name_file_R_CM= DIR+"R_CM_"+"Real"+"_"+density+"_"
		+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)
		+"_"+std::to_string(start_traj)+"_"+std::to_string(run);
		std::ofstream R_CM;
		R_CM.open(name_file_R_CM);
		
		// PRINTING INITIAL CONDITION AS CONFIGURATION 0 //

		double GR_average = 0;
		double R_EE_average = 0;
		for(std::size_t k=0;k<melt.size();k++){
			auto r_cm=compute_R_cm(melt[k]);
			R_CM<<r_cm[0]<<"\t"<<r_cm[1]<<"\t"<<r_cm[2]<<"\t";
			GR_average += gyration_radius(melt[k]);
			R_EE_average += end_end_distance(melt[k]);
			for(std::size_t z=0;z<melt[k].chain.size();z++){
				traj<<melt[k].chain[z].t<<"\t"<<melt[k].chain[z].u<<"\t"<<melt[k].chain[z].v<<"\t";
			}
		}

		GR_average = GR_average/melt.size();
		R_EE_average = R_EE_average/melt.size();
		E_bending<<E<<std::endl;
		GR<<GR_average<<std::endl;
		R_EE<<R_EE_average<<std::endl;
		traj<<std::endl;
		R_CM<<std::endl;
		traj.close();
		R_CM.close();

		long long int acc = 0;
		int N_file_to_write=MC_steps/stride;
		
		std::cout<<"CIAO"<<std::endl;

		for(int j=1;j<N_file_to_write;j++){ //loop over the N configurations to print
			for(int n_step=0;n_step<stride;n_step++){ //loop over the time between the printed configurations
				unsigned int polymer_to_move = rand()%number_of_polymers;
				unsigned int monomer_to_move = rand()%polymer_length;
				auto MC_move=rand()%12; //there are 12 possible movement
				monomer trial{melt[polymer_to_move].chain[monomer_to_move].t+mc_move[MC_move].dt,melt[polymer_to_move].chain[monomer_to_move].u+mc_move[MC_move].du,melt[polymer_to_move].chain[monomer_to_move].v+mc_move[MC_move].dv};
				if(acceptance_real_ring(lattice,trial,melt[polymer_to_move],monomer_to_move,E)){
					lattice[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]-=1;
					melt[polymer_to_move].chain[monomer_to_move].t=trial.t;
					melt[polymer_to_move].chain[monomer_to_move].u=trial.u;
					melt[polymer_to_move].chain[monomer_to_move].v=trial.v;
					melt[polymer_to_move].chain[monomer_to_move].t_image=trial.t_image;
					melt[polymer_to_move].chain[monomer_to_move].u_image=trial.u_image;
					melt[polymer_to_move].chain[monomer_to_move].v_image=trial.v_image;
					lattice[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]+=1;
					acc++;
				}
			}

			name_file3 = DIR+"traj_"+"Real"+"_"+density+"_"
			+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)
			+"_"+std::to_string(start_traj+j)+"_"+std::to_string(run);
			traj.open(name_file3);	
			
			std::string name_file_R_CM= DIR+"R_CM_"+"Real"+"_"+density+"_"
			+std::to_string(int(bending_rigidity))+"_"+std::to_string(polymer_length)+"_"
			+std::to_string(start_traj+j)+"_"+std::to_string(run);
			R_CM.open(name_file_R_CM);
			
			double GR_average = 0;
			double R_EE_average = 0;
			
			for(std::size_t k=0;k<melt.size();k++){
				auto r_cm=compute_R_cm(melt[k]);
				R_CM<<r_cm[0]<<"\t"<<r_cm[1]<<"\t"<<r_cm[2]<<"\t";
				GR_average += gyration_radius(melt[k]);
				auto R_ee_single = end_end_distance(melt[k]);
				R_EE_average += R_ee_single;
				for(std::size_t z=0;z<melt[k].chain.size();z++){
					traj<<melt[k].chain[z].t<<"\t"<<melt[k].chain[z].u<<"\t"<<melt[k].chain[z].v<<"\t";
				}
			}
			E_bending<<E<<std::endl;
			GR_average = GR_average/melt.size();
			R_EE_average = R_EE_average/melt.size();
			GR<<GR_average<<std::endl;
			R_EE<<R_EE_average<<std::endl;
			traj<<std::endl;
			R_CM<<std::endl;
			traj.close();
			R_CM.close();
		}
		
		double temp=0;
		for(auto i=0;i<Lattice_size;i++){
			for(auto j=0;j<Lattice_size;j++){
				for(auto k=0;k<Lattice_size;k++){
					if(lattice[i][j][k]>2){
						temp+=1;
						
					}
				}
			}
		}

		std::cout<<"K = "+std::to_string(int(bending_rigidity))+" ,Acceptance ratio = "<<double(acc)/double(MC_steps-stride)<<std::endl;
	}
}


#endif
