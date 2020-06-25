#pragma once
#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include <nanoflann.hpp>

#include <fstream>
#include <cstdlib>

#include "ichnos_structures.h"
#include "ichnos_utils.h"


namespace po = boost::program_options;

namespace ICHNOS {
	
	class velocityField {
	public:
		velocityField(boost::mpi::communicator& world_in);
		virtual void readVelocityField(std::string vf_file){}
		virtual void calcVelocity(vec3& vel, std::map<int, double>& proc_map, vec3& p, double& step) {}
		
		bool bIsInGhostArea(std::map<int, double> proc_map);
		int calcProcID(std::map<int, double> proc_map);
		bool InterpolateOutsideDomain = false;

	protected:
		boost::mpi::communicator world;
		double OwnerThreshold = 0.75;
		

	};

	velocityField::velocityField(boost::mpi::communicator& world_in) 
		:
		world(world_in)
	{}

	bool velocityField::bIsInGhostArea(std::map<int, double> proc_map) {
		std::map<int, double>::iterator it = proc_map.find(world.rank());
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


	//class pc3D : public velocityField {
	//public:
	//	pc3D(boost::mpi::communicator& world_in);
	//	void readVelocityField(std::string vf_file);
	//	void calcVelocity(vec3& vel, std::map<int, double>& proc_map, vec3& p);
	//	unsigned int getCloudSize() { return pcVF.pts.size(); }

	//private:
	//	pointCloud<vec3> pcVF;
	//	pointCloud<double> topElev;
	//	pointCloud<double> botElev;
	//	std::unique_ptr<nano_kd_tree_vector> velocity_tree;
	//	std::unique_ptr<nano_kd_tree_scalar> top_tree;
	//	std::unique_ptr<nano_kd_tree_scalar> bot_tree;
	//};

	//pc3D::pc3D(boost::mpi::communicator& world_in)
	//	:
	//	velocityField(world_in)
	//{}

	//void pc3D::readVelocityField(std::string vf_file) {
	//	po::options_description velocityFieldOptions("Velocity field options");
	//	po::variables_map vm_vfo;
	//	velocityFieldOptions.add_options()
	//		("pcPrefix", po::value<std::string>(), "filename for point cloud of the velocity field")
	//		("pcSuffix", po::value<std::string>(), "ending of file after procid")
	//		("OwnerThreshold", po::value<double>()->default_value(0.75), "Threshold for the processor ownership")
	//		;

	//	po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions), vm_vfo);

	//	OwnerThreshold = vm_vfo["OwnerThreshold"].as<double>();

	//	std::string prefix = vm_vfo["pcPrefix"].as<std::string>();
	//	std::string suffix;
	//	if (vm_vfo.count("pcSuffix")) {
	//		suffix = vm_vfo["pcSuffix"].as<std::string>();
	//	}
	//	else
	//		suffix = ".dat";


	//	std::string filename = prefix + std::to_string(world.rank()) + suffix;

	//	READ::readVelocityFieldFile(filename, pcVF);
	//	this->InterpolateOutsideDomain = pcVF.InterpolateOutside;
	//	this->velocity_tree = std::unique_ptr<nano_kd_tree_vector>(new nano_kd_tree_vector(3, pcVF, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
	//	velocity_tree->buildIndex();
	//
	//	/*
	//	// read and build top elevation tree
	//	std::string topfile = vm_vfo["Top"].as<std::string>();
	//	READ::read2DscalarField(topfile, topElev);
	//	this->top_tree = std::unique_ptr<nano_kd_tree_scalar>(new nano_kd_tree_scalar(3, topElev, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
	//	top_tree->buildIndex();

	//	// read and build top elevation tree
	//	std::string botfile = vm_vfo["Bot"].as<std::string>();
	//	READ::read2DscalarField(botfile, botElev);
	//	this->bot_tree = std::unique_ptr<nano_kd_tree_scalar>(new nano_kd_tree_scalar(3, botElev, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
	//	bot_tree->buildIndex();
	//	*/
	//}

	//void pc3D::calcVelocity(vec3& vel, std::map<int, double>& proc_map, vec3& p) {
	//	interpolateVectorTree(vel, proc_map, velocity_tree, p);
	//}
}
