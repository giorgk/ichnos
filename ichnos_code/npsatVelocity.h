#pragma once

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include <nanoflann.hpp>

#include "ichnos_structures.h"
#include "velocity_base.h"

namespace po = boost::program_options;
namespace ic = ICHNOS;

namespace NPSAT {
	class npsatVel : public ic::velocityField {
	public:
		npsatVel(boost::mpi::communicator& world_in);
		void readVelocityField(std::string vf_file);
		void calcVelocity(ic::vec3& vel, std::map<int, double>& proc_map, ic::vec3& p);
	private:
		ic::pointCloud<ic::vec3> VelocityCloud;
		std::unique_ptr<nano_kd_tree_vector> VelocityTree;
		ic::interpType porType;
		ic::pointCloud<double> PorosityCloud;
		double porosityValue = 1.0;
		std::unique_ptr<nano_kd_tree_scalar> PorosityTree;
		double multiplier = 1.0;
	};

	npsatVel::npsatVel(boost::mpi::communicator& world_in)
		:
		velocityField(world_in)
	{}

	void npsatVel::readVelocityField(std::string vf_file) {
		if (world.rank() == 0)
			std::cout << "Reading data..." << std::endl;
		po::options_description velocityFieldOptions("Velocity field options");
		po::variables_map vm_vfo;
		velocityFieldOptions.add_options()
			("Prefix", po::value<std::string>(), "Prefix for the filename") 
			("LeadingZeros", po::value<int>()->default_value(4), "e.g 0002->4, 000->3")
			("Suffix", po::value<std::string>(), "ending of file after procid")
			("OwnerThreshold", po::value<double>()->default_value(0.75), "Threshold for the processor ownership")
			("VelocityMultiplier", po::value<double>()->default_value(1), "This is a multiplier to scale velocity")
			("Porosity", po::value<std::string>(), "Porocity. Either a file or a single number")
			;

		po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions), vm_vfo);

		OwnerThreshold = vm_vfo["OwnerThreshold"].as<double>();
		int leadZeros = vm_vfo["LeadingZeros"].as<int>();
		multiplier = vm_vfo["VelocityMultiplier"].as<double>();

		std::string prefix = vm_vfo["Prefix"].as<std::string>();
		std::string suffix;
		if (vm_vfo.count("Suffix")) {
			suffix = vm_vfo["Suffix"].as<std::string>();
		}
		else
			suffix = ".ich";

		std::string filename = prefix + ic::num2Padstr(world.rank(), leadZeros) + suffix;
		ic::READ::readVelocityFieldFile(filename, VelocityCloud);
		this->InterpolateOutsideDomain = VelocityCloud.InterpolateOutside;
		this->VelocityTree = std::unique_ptr<nano_kd_tree_vector>(new nano_kd_tree_vector(3, VelocityCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
		VelocityTree->buildIndex();

		if (vm_vfo.count("Porosity")) {
			std::string porfile = vm_vfo["Porosity"].as<std::string>();
			if (porfile.empty()) {
				porType = ic::interpType::INGORE;
			}
			else {
				if (ic::is_input_scalar(porfile)) {
					porType = ic::interpType::SCALAR;
					porosityValue = std::stod(porfile);
				}
				else {
					porType = ic::interpType::CLOUD;
					ic::READ::read2DscalarField(porfile, PorosityCloud);
					this->PorosityTree = std::unique_ptr<nano_kd_tree_scalar>(new nano_kd_tree_scalar(3, PorosityCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
					PorosityTree->buildIndex();
				}
			}
		}
		else {
			porType = ic::interpType::INGORE;
		}
	}

	void npsatVel::calcVelocity(ic::vec3& vel, std::map<int, double>& proc_map, ic::vec3& p) {
		ic::interpolateVectorTree(vel, proc_map, VelocityTree, p);
		double porosity = 1.0;
		if (porType == ic::interpType::CLOUD)
			porosity = ic::interpolateScalarTree(PorosityTree, p);
		else if (porType == ic::interpType::SCALAR)
			porosity = porosityValue;
		vel = vel * (multiplier/porosity);
	}
}

