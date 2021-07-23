#pragma once
#include <boost/mpi.hpp>
#include <boost/program_options.hpp>

#include <fstream>
#include <cstdlib>

#include "ichnos_structures.h"
#include "ichnos_utils.h"


namespace ICHNOS {
	
	class velocityField {
	public:
		velocityField(boost::mpi::communicator& world_in);
		virtual void readVelocityField(std::string vf_file){}
		virtual void calcVelocity(vec3& vel, std::map<int, double>& proc_map, vec3& p) {}
		virtual void reset(){}
		virtual void updateStep(double& step){}
		
		bool bIsInGhostArea(std::map<int, double> proc_map);
		int calcProcID(std::map<int, double> proc_map);
		bool InterpolateOutsideDomain = true;

	protected:
		boost::mpi::communicator world;
		double OwnerThreshold = 0.75;
		

	};

	velocityField::velocityField(boost::mpi::communicator& world_in) 
		:
		world(world_in)
	{}

	bool velocityField::bIsInGhostArea(std::map<int, double> proc_map) {
		std::map<int, double>::iterator it = proc_map.find(/*dbg_rank*/world.rank());
		if (it != proc_map.end()) {
			if (it->second > OwnerThreshold)
				return false;
			else
				return true;
		}
		else {
			return true;
		}
	}

	int velocityField::calcProcID(std::map<int, double> proc_map) {
		std::map<int, double>::iterator it;
		int id = 0;
		double perc = 0;
		for (it = proc_map.begin(); it != proc_map.end(); ++it) {
			if (it->second > perc) {
				perc = it->second;
				id = it->first;
			}
		}
		return id;
	}
}
