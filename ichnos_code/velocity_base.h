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
		virtual bool readVelocityField(std::string vf_file, int nPnts){return true;}
		virtual void calcVelocity(vec3& vel,
                                  std::vector<int>& ids,
                                  std::vector<double>& weights,
                                  double tm = 0) {}
		virtual void reset(){}
		virtual void updateStep(double& step){}
		virtual void getVec3Data(std::vector<vec3>& data){}
		
		bool bIsInGhostArea(std::map<int, double> proc_map);
		VelType getVelType(){return Vtype;};
		int calcProcID(std::map<int, double> proc_map);
		bool InterpolateOutsideDomain = true;
		void SetStepOptions(StepOptions stepOpt_in);

	protected:
		boost::mpi::communicator world;
		double OwnerThreshold = 0.75;
		VelType Vtype;
		StepOptions stepOpt;
		

	};

	velocityField::velocityField(boost::mpi::communicator& world_in)
		:
		world(world_in)
	{}

	void velocityField::SetStepOptions(StepOptions stepOpt_in) {
        stepOpt = stepOpt_in;
	}

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
