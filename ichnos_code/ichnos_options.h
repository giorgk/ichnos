#pragma once

#include <iostream>

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>

#include "ichnos_structures.h"


namespace po = boost::program_options;


namespace ICHNOS {
	enum class VelType {
	//	Cloud3d,
		IWFM,
		NPSAT,
		STOCH,
		INVALID
	};

	std::string castVelType2String(VelType vt) {
		std::map <VelType, std::string> vtMap;
		std::map <VelType, std::string>::iterator it;
		//vtMap.insert(std::pair<VelType, std::string>(VelType::Cloud3d, "Cloud3d"));
		vtMap.insert(std::pair<VelType, std::string>(VelType::IWFM, "IWFM"));
		vtMap.insert(std::pair<VelType, std::string>(VelType::NPSAT, "NPSAT"));
		vtMap.insert(std::pair<VelType, std::string>(VelType::STOCH, "STOCH"));
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
		//vtMap.insert(std::pair<std::string, VelType>("Cloud3d", VelType::Cloud3d));
		vtMap.insert(std::pair<std::string, VelType>("IWFM", VelType::IWFM));
		vtMap.insert(std::pair<std::string, VelType>("NPSAT", VelType::NPSAT));
		vtMap.insert(std::pair<std::string, VelType>("STOCH", VelType::STOCH));
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
		int nproc = 0;
		int niter = 0;

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
				std::cout << "| Version : 0.1.10 |" << std::endl;
				std::cout << "|    by  giorgk    |" << std::endl;
				std::cout << "|------------------|" << std::endl;
			}
			return false;
		}

		// Configuration file options
		po::options_description config_options("Configuration file options");
		config_options.add_options()
			// Velocity Options
			("Velocity.Type", po::value<std::string>(), "Type of velocity. (NPSAT or IWFM)")
			("Velocity.ConfigFile", po::value<std::string >(), "Set configuration file for the velocity field")

			// Domain options
			("Domain.Outline", po::value<std::string >(), "A filename that containts the vertices of the outline polygon")
			("Domain.TopFile", po::value<std::string >(), "A filename with the point cloud of the top elevation")
			("Domain.TopRadius", po::value<double>()->default_value(1000), "Search Radious for top elevation")
			("Domain.TopPower", po::value<double>()->default_value(3), "Search Power for top elevation")
			("Domain.BottomFile", po::value<std::string >(), "A filename with the point cloud of the bottom elevation")
			("Domain.BottomRadius", po::value<double>()->default_value(1000), "Search Radious for bottom elevation")
			("Domain.BottomPower", po::value<double>()->default_value(3), "Search Power for bottom elevation")
			("Domain.ProcessorPolys", po::value<std::string >(), "A filename that containts the coordinates of each processor polygon")
			("Domain.ExpandedPolys", po::value<std::string >(), "A filename that containts the coordinates of each Expanded polygon")
			("Domain.AttractFile", po::value<std::string >(), "A filename that containts the particle attractors")
			("Domain.AttractRadius", po::value<double>()->default_value(50), "Particles closer to this distance will be captured from attractors")

			// Stopping criteria
			("StoppingCriteria.MaxIterationsPerStreamline", po::value<int>()->default_value(1000), "Maximum number of steps per streamline")
			("StoppingCriteria.MaxProcessorExchanges", po::value<int>()->default_value(50), "Maximum number of that a particles are allowed to change processors")
			("StoppingCriteria.AgeLimit", po::value<double>()->default_value(-1), "If the particle exceed this limit stop particle tracking. Negative means no limit")
			("StoppingCriteria.Stuckiter", po::value<int>()->default_value(10), "After Stuckiter exit particle tracking")

			// Step configuration
			("StepConfig.Method", po::value<std::string >(), "Method for steping")
			("StepConfig.Direction", po::value<double>()->default_value(1), "Backward or forward particle tracking")
			("StepConfig.StepSize", po::value<double>()->default_value(1), "Step Size in units of length")
			("StepConfig.minExitStepSize", po::value<double>()->default_value(0.1), "Minimum Step Size at the exit as percentage of the stepsize")

			// Adaptive Step configurations
			("AdaptStep.UpdateStepSize", po::value<int>()->default_value(1), "Update step size wrt bbox")
			("AdaptStep.MaxStepSize", po::value<double>()->default_value(2), "Maximum Step Size in units of length")
			("AdaptStep.MinStepSize", po::value<double>()->default_value(0.1), "Minimum Step Size in units of length")
			("AdaptStep.increasRatechange", po::value<double>()->default_value(1.5), "Maximum Step Size in units of length")
			("AdaptStep.limitUpperDecreaseStep", po::value<double>()->default_value(0.75), "Upper limit of decrease step size")
			("AdaptStep.Tolerance", po::value<double>()->default_value(0.1), "Tolerance when the RK45 is used")

			// InputOutput
			("InputOutput.PartilceFile", po::value<std::string >(), "A filename with the initial positions of particles")
			("InputOutput.WellFile", po::value<std::string >(), "A filename with the well locations")
			("InputOutput.OutputFile", po::value<std::string >(), "Prefix for the output file")
			("InputOutput.ParticlesInParallel", po::value<int>()->default_value(1000), "Maximum number run in parallel")
			("InputOutput.GatherOneFile", po::value<int>()->default_value(1), "Put all streamlines into one file")

			// Other
			("Other.Nrealizations", po::value<int>()->default_value(1), "NUmber of realizations")

			// Unsued
			("Unused.nThreads", po::value<int>()->default_value(1), "Number of threads")
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
			try {
				Popt.configfile = vm_cmd["config"].as<std::string>().c_str();
				std::cout << "Configuration file: " << vm_cmd["config"].as<std::string>().c_str() << std::endl;
				po::store(po::parse_config_file<char>(vm_cmd["config"].as<std::string>().c_str(), config_options), vm_cfg);

				{// Velocity Options
					velocityFieldType = castVelType2Enum(vm_cfg["Velocity.Type"].as<std::string>());
					if (velocityFieldType == VelType::INVALID) {
						std::cout << vm_cfg["Velocity.Type"].as<std::string>() << " is an invalid velocity type" << std::endl;
						return false;
					}
					velocityFieldFileName = vm_cfg["Velocity.ConfigFile"].as<std::string>();
				}

				{// Domain options
					Dopt.polygonFile = vm_cfg["Domain.Outline"].as<std::string>();
					// Top elevation parameters
					Dopt.TopElevationFile = vm_cfg["Domain.TopFile"].as<std::string>();
					Dopt.TopRadius = vm_cfg["Domain.TopRadius"].as<double>();
					Dopt.TopRadius = Dopt.TopRadius*Dopt.TopRadius;
					Dopt.TopPower = vm_cfg["Domain.TopPower"].as<double>();
					// Bottom elevation parameters
					Dopt.BottomeElevationFile = vm_cfg["Domain.BottomFile"].as<std::string>();
					Dopt.BotRadius = vm_cfg["Domain.BottomRadius"].as<double>();
					Dopt.BotRadius = Dopt.BotRadius*Dopt.BotRadius;
					Dopt.BotPower = vm_cfg["Domain.BottomPower"].as<double>();
					Dopt.processorDomainFile = vm_cfg["Domain.ProcessorPolys"].as<std::string>();
					Dopt.expandedDomainFile = vm_cfg["Domain.ExpandedPolys"].as<std::string>();
					Dopt.myRank = world.rank();
					Dopt.AttractorsFile = vm_cfg["Domain.AttractFile"].as<std::string>();
					Dopt.AttractRadius = vm_cfg["Domain.AttractRadius"].as<double>();
					Dopt.AttractRadius = Dopt.AttractRadius * Dopt.AttractRadius;
				}

				{// Stopping criteria
					Popt.MaxIterationsPerStreamline = vm_cfg["StoppingCriteria.MaxIterationsPerStreamline"].as<int>();
					Popt.MaxProcessorExchanges = vm_cfg["StoppingCriteria.MaxProcessorExchanges"].as<int>();
					Popt.AgeLimit = vm_cfg["StoppingCriteria.AgeLimit"].as<double>();
					Popt.StuckIterations = vm_cfg["StoppingCriteria.Stuckiter"].as<int>();
				}

				{// Step configuration
					SolutionMethods method = castMethod2Enum(vm_cfg["StepConfig.Method"].as<std::string>());
					if (method == SolutionMethods::INVALID)
						return false;
					else {
						Popt.method = method;
					}
					double tmp = vm_cfg["StepConfig.Direction"].as<double>();
					if (tmp >= 0)
						Popt.Direction = 1;
					else
						Popt.Direction = -1;
					
					Popt.StepSize = vm_cfg["StepConfig.StepSize"].as<double>();

					Popt.minExitStepSize = vm_cfg["StepConfig.minExitStepSize"].as<double>();
					if (Popt.minExitStepSize < 0 || Popt.minExitStepSize > 1) {
						std::cout << "minExitStepSize should be between 0 and 1. It gets the default value of 0.1" << std::endl;
						Popt.minExitStepSize = 0.1;
					}
				}

				{// Adaptive Step configurations
					Popt.UpdateStepSize = vm_cfg["AdaptStep.UpdateStepSize"].as<int>();
					Popt.MaxStepSize = vm_cfg["AdaptStep.MaxStepSize"].as<double>();
					Popt.MinStepSize = vm_cfg["AdaptStep.MinStepSize"].as<double>();
					Popt.increasRatechange = vm_cfg["AdaptStep.increasRatechange"].as<double>();
					if (Popt.increasRatechange < 1) {
						std::cout << "increasRatechange should be higher than 1. It gets the default value of 1.5" << std::endl;
						Popt.increasRatechange = 1.5;
					}
					Popt.limitUpperDecreaseStep = vm_cfg["AdaptStep.limitUpperDecreaseStep"].as<double>();
					if (Popt.limitUpperDecreaseStep < 0 || Popt.limitUpperDecreaseStep > 1) {
						std::cout << "limitUpperDecreaseStep should be between 0 and 1. It gets the default value of 0.75" << std::endl;
						Popt.limitUpperDecreaseStep = 0.75;
					}
					Popt.ToleranceStepSize = vm_cfg["AdaptStep.Tolerance"].as<double>();
				}

				{// InputOutput
					Popt.ParticleFile = vm_cfg["InputOutput.PartilceFile"].as<std::string>();
					Popt.WellFile = vm_cfg["InputOutput.WellFile"].as<std::string>();
					Popt.OutputFile = vm_cfg["InputOutput.OutputFile"].as<std::string>();
					Popt.ParticlesInParallel = vm_cfg["InputOutput.ParticlesInParallel"].as<int>(); 
				}

				{// Other
					Popt.Nrealizations = vm_cfg["Other.Nrealizations"].as<int>();
				}

				{// Unsued
					nThreads = vm_cfg["Unused.nThreads"].as<int>();
				}
			
			}
			catch (std::exception& E)
			{
				std::cout << E.what() << std::endl;
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
				std::cout << "In gather mode you have to specify the number of iterations (iter) used in particle tracking" << std::endl;
				return false;
			}
			{
				int tmp = vm_cfg["InputOutput.GatherOneFile"].as<int>();
				if (tmp == 0)
					Popt.gatherOneFile = false;
				else
					Popt.gatherOneFile = true;
			}
		}
		return true;
	}
}
