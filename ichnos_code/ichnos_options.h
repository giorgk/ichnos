#pragma once

#include <iostream>

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>

#include "ichnos_structures.h"


namespace po = boost::program_options;


namespace ICHNOS {
	enum class VelType {
		Cloud3d,
		IWFM,
		INVALID
	};

	std::string castVelType2String(VelType vt) {
		std::map <VelType, std::string> vtMap;
		std::map <VelType, std::string>::iterator it;
		vtMap.insert(std::pair<VelType, std::string>(VelType::Cloud3d, "Cloud3d"));
		vtMap.insert(std::pair<VelType, std::string>(VelType::IWFM, "IWFM"));
		it = vtMap.find(vt);
		if (it != vtMap.end())
			return it->second;
		else {
			return "INVALID";
		}
	}

	VelType castVelType2Enum(std::string vt) {
		std::map < std::string, VelType> vtMap;
		std::map < std::string, VelType>::iterator it;
		vtMap.insert(std::pair<std::string, VelType>("Cloud3d", VelType::Cloud3d));
		vtMap.insert(std::pair<std::string, VelType>("IWFM", VelType::IWFM));
		it = vtMap.find(vt);
		if (it != vtMap.end())
			return it->second;
		else {
			return VelType::INVALID;
		}
	}

	class options {
	public:
		options(boost::mpi::communicator& world_in);

		bool readInput(int argc, char* argv[]);

		std::string getVelFname() {
			return velocityFieldFileName;
		}

		ParticleOptions Popt;
		DomainOptions Dopt;
		bool gatherMode = false;
		int nproc;
		int niter;

		VelType velocityFieldType;
		

		

		
	private:
		boost::mpi::communicator world;
		bool bIsMultiThreaded = true;
		int nThreads = 1;
		std::string velocityFieldFileName;
	};

	options::options(boost::mpi::communicator& world_in)
		:
		world(world_in)
	{
		if (world.size() > 1) {
			bIsMultiThreaded = false;
		}
	}

	bool options::readInput(int argc, char* argv[]) {
		// Command line options
		po::options_description commandLineOptions("Command line options");
		commandLineOptions.add_options()
			("version,v", "print version information")
			("help,h", "Get a list of options in the configuration file")
			("config,c", po::value<std::string >(), "Set configuration file")
			("gather,g", "Triggers the gather mode. The configuration file is still needed")
			("nproc,n", po::value<int>(), "number of processors (used in gather mode only)")
			("iter,i", po::value<int>(), "number of iterations (used in gather mode only)")
			;

		po::variables_map vm_cmd;

		po::store(po::parse_command_line(argc, argv, commandLineOptions), vm_cmd);

		if (vm_cmd.size() == 0) {
			if (world.rank() == 0) {
				std::cout << " To run Ichnos specify the configuration file as" << std::endl;
				std::cout << "-c config" << std::endl << std::endl;;
				std::cout << "Other command line options are:" << std::endl;
				std::cout << commandLineOptions << std::endl;
			}
			return false;
		}

		if (vm_cmd.count("version")) {
			if (world.rank() == 0) {
				std::cout << "|------------------|" << std::endl;
				std::cout << "|      ICHNOS      |" << std::endl;
				std::cout << "| Version : 0.0.02 |" << std::endl;
				std::cout << "|    by  giorgk    |" << std::endl;
				std::cout << "|------------------|" << std::endl;
			}
			return false;
		}

		// Configuration file options
		po::options_description config_options("Configuration file options");
		config_options.add_options()
			("nThreads", po::value<int>()->default_value(1), "Number of threads")
			("VelocityConfig", po::value<std::string >(), "Set configuration file for the velocity field")
			("Stuckiter", po::value<int>()->default_value(10), "After Stuckiter exit particle tracking")
			("Method", po::value<std::string >(), "Method for time steping")
			("DomainPolygon", po::value<std::string >(), "A filename that containts the vertices of the outline polygon")
			("TopElevation", po::value<std::string >(), "A filename with the point cloud of the top elevation")
			("BottomElevation", po::value<std::string >(), "A filename with the point cloud of the bottom elevation")
			("PartilceFile", po::value<std::string >(), "A filename with the initial positions of particles")
			("WellFile", po::value<std::string >(), "A filename with the well locations")
			("Direction", po::value<double>()->default_value(1), "Backward or forward particle tracking")
			("StepSize", po::value<double>()->default_value(1), "Step Size in units of length")
			("MaxStepSize", po::value<double>()->default_value(2), "Maximum Step Size in units of length")
			("increasRatechange", po::value<double>()->default_value(1.5), "Maximum Step Size in units of length")
			("MinStepSize", po::value<double>()->default_value(0.1), "Minimum Step Size in units of length")
			("minExitStepSize", po::value<double>()->default_value(0.1), "Minimum Step Size at the exit as percentage of the stepsize")
			("limitUpperDecreaseStep", po::value<double>()->default_value(0.75), "Upper limit of decrease step size")
			("ToleranceStepSize", po::value<double>()->default_value(0.1), "Tolerance when the RK45 is used")
			("MaxIterationsPerStreamline", po::value<int>()->default_value(1000), "Maximum number of steps per streamline")
			("MaxProcessorExchanges", po::value<int>()->default_value(50), "Maximum number of that a particles are allowed to change processors")
			("OutputFile", po::value<std::string >(), "Prefix for the output file")
			("ParticlesInParallel", po::value<int>()->default_value(1000), "Maximum number run in parallel")
			("GatherOneFile", po::value<int>()->default_value(1), "Put all streamlines into one file")
			("VelocityType", po::value<std::string>(), "Type of velocity. (Cloud3d or IWFM)")
			;

		if (vm_cmd.count("help")) {
			if (world.rank() == 0) {
				std::cout << " To run ICHNOS specify the configuration file as" << std::endl;
				std::cout << "-c config" << std::endl << std::endl;;
				std::cout << "Other command line options are:" << std::endl;
				std::cout << commandLineOptions << std::endl;

				std::cout << "ICHNOS configuration file options:" << std::endl;
				std::cout << "The options without default values are mandatory" << std::endl;
				std::cout << "(All options are case sensitive)" << std::endl;
				std::cout << "------------------------------" << std::endl;
				std::cout << config_options << std::endl;
			}
			return false;
		}

		po::variables_map vm_cfg;
		if (vm_cmd.count("config")) {
			Popt.configfile = vm_cmd["config"].as<std::string>().c_str();
			std::cout << "Configuration file: " << vm_cmd["config"].as<std::string>().c_str() << std::endl;
			po::store(po::parse_config_file<char>(vm_cmd["config"].as<std::string>().c_str(), config_options), vm_cfg);
			nThreads = vm_cfg["nThreads"].as<int>();
			velocityFieldFileName = vm_cfg["VelocityConfig"].as<std::string>();

			// particle tracking options
			Popt.StuckIterations = vm_cfg["Stuckiter"].as<int>();
			SolutionMethods method = castMethod2Enum(vm_cfg["Method"].as<std::string>());
			if (method == SolutionMethods::INVALID)
				return false;
			else {
				Popt.method = method;
			}
			double tmp = vm_cfg["Direction"].as<double>();
			if (tmp >= 0)
				Popt.Direction = 1;
			else
				Popt.Direction = -1;

			Popt.StepSize = vm_cfg["StepSize"].as<double>();
			Popt.MaxStepSize = vm_cfg["MaxStepSize"].as<double>();
			Popt.MinStepSize = vm_cfg["MinStepSize"].as<double>();
			Popt.increasRatechange = vm_cfg["increasRatechange"].as<double>();
			if (Popt.increasRatechange < 1) {
				std::cout << "increasRatechange should be higher than 1. It gets the default value of 1.5" << std::endl;
				Popt.increasRatechange = 1.5;
			}
			Popt.limitUpperDecreaseStep = vm_cfg["limitUpperDecreaseStep"].as<double>();
			if (Popt.limitUpperDecreaseStep < 0 || Popt.limitUpperDecreaseStep > 1) {
				std::cout << "limitUpperDecreaseStep should be between 0 and 1. It gets the default value of 0.75" << std::endl;
				Popt.limitUpperDecreaseStep = 0.75;
			}
			Popt.minExitStepSize = vm_cfg["minExitStepSize"].as<double>();
			if (Popt.minExitStepSize < 0 || Popt.minExitStepSize > 1) {
				if (Popt.minExitStepSize < 0 || Popt.minExitStepSize > 1) {
					std::cout << "minExitStepSize should be between 0 and 1. It gets the default value of 0.1" << std::endl;
					Popt.minExitStepSize = 0.1;
				}
			}
			Popt.minExitStepSize = Popt.StepSize * Popt.minExitStepSize;
			Popt.ToleranceStepSize = vm_cfg["ToleranceStepSize"].as<double>();
			Popt.MaxIterationsPerStreamline = vm_cfg["MaxIterationsPerStreamline"].as<int>();
			Popt.MaxProcessorExchanges = vm_cfg["MaxProcessorExchanges"].as<int>();
			Popt.ParticlesInParallel = vm_cfg["ParticlesInParallel"].as<int>(); 
			Popt.ParticleFile = vm_cfg["PartilceFile"].as<std::string>();
			Popt.WellFile = vm_cfg["WellFile"].as<std::string>();
			Popt.OutputFile = vm_cfg["OutputFile"].as<std::string>();

			// Domain options
			Dopt.polygonFile = vm_cfg["DomainPolygon"].as<std::string>();
			Dopt.TopElevationFile = vm_cfg["TopElevation"].as<std::string>();
			Dopt.BottomeElevationFile = vm_cfg["BottomElevation"].as<std::string>();

			velocityFieldType = castVelType2Enum(vm_cfg["VelocityType"].as<std::string>());
			if (velocityFieldType == VelType::INVALID) {
				return false;
			}
		}

		

		if (vm_cmd.count("gather")) {
			gatherMode = true;
			if (vm_cmd.count("nproc"))
				nproc = vm_cmd["nproc"].as<int>();
			else {
				std::cout << "In gather mode you have to specify the number of processors (nproc) used in particle tracking" << std::endl;
				return false;
			}
			if (vm_cmd.count("iter"))
				niter = vm_cmd["iter"].as<int>();
			else {
				std::cout << "In gather mode you have to specify the number of iterations (iter) used particle tracking" << std::endl;
				return false;
			}
			{
				int tmp = vm_cfg["GatherOneFile"].as<int>();
				if (tmp == 0)
					Popt.gatherOneFile = false;
				else
					Popt.gatherOneFile = true;
			}
		}
		return true;
	}
}
