#pragma once
#include <boost/mpi.hpp>
#include <boost/program_options.hpp>

#include <fstream>
#include <cstdlib>

#include "ichnos_structures.h"
#include "ichnos_utils.h"
#include "ichnos_XYZ_base.h"


namespace ICHNOS {
	
	class velocityField {
	public:
		velocityField(boost::mpi::communicator& world_in, XYZ_base &XYZ_in);
		virtual bool readVelocityField(std::string vf_file){return true;}
		virtual void calcVelocity(vec3& p, vec3& vel,
                                  std::map<int, double>& proc_map,
                                  helpVars& pvlu,
                                  bool& out,
                                  //std::vector<int> &ids,
                                  //std::vector<double> &weights,
                                  double tm = 0) {}
		virtual void reset(Streamline& S){}
		virtual double stepTimeupdate(helpVars& pvlu){return 0;}
		virtual void updateStep(double& step){}
		virtual void getVec3Data(std::vector<vec3>& data){}
		
		bool bIsInGhostArea(std::map<int, double> proc_map);
		VelType getVelType(){return Vtype;}
        bool isVelTransient(){return isVeltrans;}
		int calcProcID(std::map<int, double> proc_map);
		bool InterpolateOutsideDomain = true;
		void SetStepOptions(StepOptions stepOpt_in);
        XYZ_base& XYZ;

	protected:
		boost::mpi::communicator world;
		double OwnerThreshold = 0.75;
		VelType Vtype;
		StepOptions stepOpt;
		std::string Prefix;
		std::string Suffix;
		int leadingZeros;
        double multiplier = 1.0;
        bool isVeltrans;

	};

	velocityField::velocityField(boost::mpi::communicator& world_in, XYZ_base &XYZ_in)
		:
		world(world_in),
        XYZ(XYZ_in)
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
