#ifndef __initial_condition__
#define __initial_condition__

#include <fstream>
#include <string>
#include <iostream>
#include "polymer.hpp"
#include "initialize_MC_moves.hpp"

extern std::map<std::vector<int>,int> MC_move_n;

void Initial_Condition(int polymer_length,int number_of_polymers){
		
		double monomer_density= double(polymer_length*number_of_polymers)/double(Lattice_size*Lattice_size*Lattice_size/2);
		std::cout<<"Density = "<< monomer_density <<std::endl;
		std::string density;
		if(number_of_polymers==1){
			density="0";
		}
		else{
			auto d = int(monomer_density*100); 
			density= std::to_string(d);
		}
		
		int lattice[Lattice_size][Lattice_size][Lattice_size]={0};
		
		std::vector<MC_move> mc_move=initialize_MC_move();

		int number_of_created_polymer=0;
		
		std::vector<polymer> melt;
		
		int stored_length=1;
		
		while(number_of_created_polymer < number_of_polymers){
			polymer polymer_tmp{polymer_length,number_of_created_polymer};
			
			int t= rand() % Lattice_size;
			int u= rand() % Lattice_size;
			int v= rand() % Lattice_size;
			
			while(lattice[t][u][v] != 0 || (t+u+v)%2 != 0){
				t= rand() % Lattice_size;
				u= rand() % Lattice_size;
				v= rand() % Lattice_size;
			}

			//std::cout<<number_of_created_polymer<<"\t"<<lattice[t][u][v]<<"\t"<<(t+u+v)%2<<"\t"<<t<<"\t"<<u<<"\t"<<v<<std::endl;
			
			monomer monomer_tmp{t,u,v};
			polymer_tmp.chain.push_back(monomer_tmp);
			polymer_tmp.chain.push_back(monomer_tmp);
			lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]+=2;
			int n_monomers=2;
			int occupation=0;

			for(int j=1;j<int(polymer_length);j++){
				int n_trials=1;
				auto MC_move = rand()%12; //there are 12 possible movement
				monomer monomer_tmp{t+mc_move[MC_move].dt,u+mc_move[MC_move].du,v+mc_move[MC_move].dv};
				occupation = lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image];
				while(occupation != 0 && n_trials < 30){
					auto MC_move = rand()%12; //there are 12 possible movement
					monomer monomer_trial{t+mc_move[MC_move].dt,u+mc_move[MC_move].du,v+mc_move[MC_move].dv};
					monomer_tmp = monomer_trial;
					occupation = lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image];
					n_trials+=1;
				}

				if(n_trials <= 30 && occupation==0){
					t=monomer_tmp.t;
					u=monomer_tmp.u;
					v=monomer_tmp.v;
					polymer_tmp.chain.push_back(monomer_tmp);
					polymer_tmp.chain.push_back(monomer_tmp);
					lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]+=2;
					n_monomers+=2;
				}

				else{
					//std::cout<<polymer_tmp.chain.size()<<std::endl;
					for(std::size_t i=0;i<polymer_tmp.chain.size();i++){
						lattice[polymer_tmp.chain[i].t_image][polymer_tmp.chain[i].u_image][polymer_tmp.chain[i].v_image]-=2;
					}
					//std::cout<<j<<std::endl;
					break;
				}
			
			if(n_monomers==int(polymer_length)){
				melt.push_back(polymer_tmp);
				number_of_created_polymer+=1;
				std::cout<<number_of_created_polymer<<std::endl;
			}

			}
		}

		std::cout<<"END"<<std::endl;

		
		bool condition=true;

		for(auto i=0;i<Lattice_size;i++){
			for(auto j=0;j<Lattice_size;j++){
				for(auto k=0;k<Lattice_size;k++){
					if(lattice[i][j][k] !=1 && lattice[i][j][k] !=0){
						condition=false;
						std::cout<<lattice[i][j][k]<<std::endl;
						std::cout<<i<<"\t"<<j<<"\t"<<k<<std::endl;
					}
				}
			}
		}
		
		for(std::size_t k=0;k<melt.size();k++){
				for(std::size_t z=0;z<melt[k].chain.size()-1;z++){
					if(fabs(distance(melt[k].chain[z],melt[k].chain[z+1])-1)>0.1){
						condition=false;	
						std::cout<<fabs(distance(melt[k].chain[z],melt[k].chain[z+1]))<<"\t"<<k<<"\t"<<z<<std::endl;
					}
				}
		}
		
		if(condition){
			std::string name_file;
			if(density=="0"){
				name_file = "../ICS/IC_"+std::to_string(polymer_length);
			}
			else{
				name_file = "../ICS/IC_"+density+"_"+std::to_string(polymer_length);
			}

			std::cout<<name_file<<std::endl;
			std::ofstream traj;
			traj.open(name_file);

			for(std::size_t k=0;k<melt.size();k++){
				for(std::size_t z=0;z<melt[k].chain.size();z++){
					traj<<melt[k].chain[z].t<<"\t"<<melt[k].chain[z].u<<"\t"<<melt[k].chain[z].v<<"\t";
				}
			}
			traj.close();
		}

		else{
			std::cout<<"There's been an error somewhere during the creation of the initial condition"<<std::endl;
			std::string name_file = "IC_"+density+"_"+std::to_string(polymer_length)+"_error";
			std::ofstream traj;
			traj.open(name_file);
			traj<<"There's been an error somewhere during the creation of the initial condition"<<std::endl;
		}
}

void Initial_Condition_Ring(int polymer_length,int number_of_polymers){

		std::cout<<polymer_length<<"\t"<<number_of_polymers<<std::endl;

		double monomer_density= double(polymer_length*number_of_polymers)/double(Lattice_size*Lattice_size*Lattice_size/2);
		std::cout<<monomer_density<<std::endl;
		std::string length= std::to_string(int(monomer_density*100));

		int lattice[Lattice_size][Lattice_size][Lattice_size]={0};
		
		std::vector<MC_move> mc_move=initialize_MC_move();
		
		int number_of_created_polymer=0;
		
		std::vector<polymer> melt;

		while(number_of_created_polymer < number_of_polymers){
			
			polymer polymer_tmp{polymer_length,number_of_created_polymer};
			
			int t= rand() % Lattice_size;
			int u= rand() % Lattice_size;
			int v= rand() % Lattice_size;
			
			while(lattice[t][u][v] != 0 || (t+u+v)%2 != 0){
				t= rand() % Lattice_size;
				u= rand() % Lattice_size;
				v= rand() % Lattice_size;
			}

			std::cout<<number_of_created_polymer<<"\t"<<lattice[t][u][v]<<"\t"<<(t+u+v)%2<<"\t"<<t<<"\t"<<u<<"\t"<<v<<std::endl;
			
			monomer monomer_tmp{t,u,v};
			polymer_tmp.chain.push_back(monomer_tmp);
			lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]+=1;
			int n_monomers=1;
			int occupation=0;

			for(int j=1;j<int(polymer_length/2);j++){
				
				int n_trials=1;
				auto MC_move = rand()%12; //there are 12 possible movement
				monomer monomer_tmp{t+mc_move[MC_move].dt,u+mc_move[MC_move].du,v+mc_move[MC_move].dv};
				occupation = lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image];

				while(occupation != 0 && n_trials < 12){
					auto MC_move = rand()%12; //there are 12 possible movement
					monomer monomer_trial{t+mc_move[MC_move].dt,u+mc_move[MC_move].du,v+mc_move[MC_move].dv};
					monomer_tmp = monomer_trial;
					occupation = lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image];
					n_trials+=1;
				}

				if(n_trials <= 12 && occupation ==0){
					t=monomer_tmp.t;
					u=monomer_tmp.u;
					v=monomer_tmp.v;
					polymer_tmp.chain.push_back(monomer_tmp);
					lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]+=1;
					n_monomers+=1;
				}

				else{
					
					for(std::size_t i=0;i<polymer_tmp.chain.size();i++){
						lattice[polymer_tmp.chain[i].t_image][polymer_tmp.chain[i].u_image][polymer_tmp.chain[i].v_image]-=1;
					}
					//std::cout<<"break at n = "<<number_of_created_polymer<<std::endl;
					break;
				}
			
			}

			if(n_monomers==int(polymer_length/2)){
				for(int i=n_monomers-1;i>=0;i--){
					polymer_tmp.chain.push_back(polymer_tmp.chain[i]);
					n_monomers+=1;
					lattice[polymer_tmp.chain[i].t_image][polymer_tmp.chain[i].u_image][polymer_tmp.chain[i].v_image]+=1;
				}
				melt.push_back(polymer_tmp);
				number_of_created_polymer+=1;
			}
		}
		
		bool condition=true;

		for(auto i=0;i<Lattice_size;i++){
			for(auto j=0;j<Lattice_size;j++){
				for(auto k=0;k<Lattice_size;k++){
					if(lattice[i][j][k] != 2 && lattice[i][j][k] != 0){
						condition=false;
						std::cout<<lattice[i][j][k]<<std::endl;
						std::cout<<i<<"\t"<<j<<"\t"<<k<<std::endl;
					}
				}
			}
		}

		
		if(condition){
			for(std::size_t k=0;k<melt.size();k++){
				for(std::size_t z=0;z<melt[k].chain.size()-1;z++){
					if(fabs(distance(melt[k].chain[z],melt[k].chain[z+1])-1)>0.1 && z != int(melt[k].chain.size()/2)-1){
						condition=false;
						std::cout<<distance(melt[k].chain[z],melt[k].chain[z+1])<<std::endl;
					}
				}
			}
		}
		
		if(condition){
			std::string name_file = "../ICS_RINGS/IC_RING_"+length+"_"+std::to_string(polymer_length);
			std::cout<<name_file<<std::endl;
			std::ofstream traj;
			traj.open(name_file);

			for(std::size_t k=0;k<melt.size();k++){
				for(std::size_t z=0;z<melt[k].chain.size();z++){
					traj<<melt[k].chain[z].t<<"\t"<<melt[k].chain[z].u<<"\t"<<melt[k].chain[z].v<<"\t";
				}
			}
			traj.close();
		}

		else{
			std::cout<<"There's been an error somewhere during the creation of the initial condition"<<std::endl;
			std::string name_file = "IC_"+length+"_"+std::to_string(polymer_length)+"_error";
			std::ofstream traj;
			traj.open(name_file);
			traj<<"There's been an error somewhere during the creation of the initial condition"<<std::endl;
		}
}

#endif

