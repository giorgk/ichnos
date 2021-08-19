#pragma once

#include <iostream>

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>

#include "ichnos_structures.h"


namespace po = boost::program_options;


namespace ICHNOS {
    /*!
     * VelType is a list of velocity field types
     */
	enum class VelType {
		STEADY, /// This is a velocity type where a single value for the velocity is defined
		TRANS, /// This is a velocity type where for its point the velocity is defined as a time series
		STOCH, /// (Experimental) This is a velocity type where the velocity is defined in some stochastic manner
		INVALID /// This is used for any velocity type that does not make any sense
	};

	/*!
	 * A utility function that converts the enumeration to string type
	 * @param vt The enumeration to convert
	 * @return The string representation of the enumeration
	 */
	std::string castVelType2String(VelType vt) {
		std::map <VelType, std::string> vtMap;
		std::map <VelType, std::string>::iterator it;
		//vtMap.insert(std::pair<VelType, std::string>(VelType::Cloud3d, "Cloud3d"));
        vtMap.insert(std::pair<VelType, std::string>(VelType::STEADY, "STEADY"));
		vtMap.insert(std::pair<VelType, std::string>(VelType::TRANS, "TRANS"));
		vtMap.insert(std::pair<VelType, std::string>(VelType::STOCH, "STOCH"));
		it = vtMap.find(vt);
		if (it != vtMap.end())
			return it->second;
		else {
			return "INVALID";
		}
	}
    /*!
     * A utility function that converts a string to VelType enumeration
     * @param vt the string to convert
     * @return The corresponing enumeration. Retunds VelType::INVALID if the input is not recognized
     */
	VelType castVelType2Enum(std::string vt) {
		std::map < std::string, VelType> vtMap;
		std::map < std::string, VelType>::iterator it;
		//vtMap.insert(std::pair<std::string, VelType>("Cloud3d", VelType::Cloud3d));
		vtMap.insert(std::pair<std::string, VelType>("STEADY", VelType::STEADY));
        vtMap.insert(std::pair<std::string, VelType>("TRANS", VelType::TRANS));
		vtMap.insert(std::pair<std::string, VelType>("STOCH", VelType::STOCH));
		it = vtMap.find(vt);
		if (it != vtMap.end())
			return it->second;
		else {
			return VelType::INVALID;
		}
	}

	enum class XYZType{
	    CLOUD,
	    IWFM,
	    INVALID
	};

	std::string castXYZType2String(XYZType xyzt) {
	    std::map <XYZType, std::string> xyztMap;
	    std::map <XYZType, std::string>::iterator it;
	    //vtMap.insert(std::pair<VelType, std::string>(VelType::Cloud3d, "Cloud3d"));
	    xyztMap.insert(std::pair<XYZType, std::string>(XYZType::CLOUD, "CLOUD"));
	    xyztMap.insert(std::pair<XYZType, std::string>(XYZType::IWFM, "IWFM"));
	    it = xyztMap.find(xyzt);
	    if (it != xyztMap.end())
	        return it->second;
	    else {
	        return "INVALID";
	    }
	}

	XYZType castXYZType2Enum(std::string xyzt) {
	    std::map < std::string, XYZType> xyztMap;
	    std::map < std::string, XYZType>::iterator it;
	    //vtMap.insert(std::pair<std::string, VelType>("Cloud3d", VelType::Cloud3d));
	    xyztMap.insert(std::pair<std::string, XYZType>("CLOUD", XYZType::CLOUD));
	    xyztMap.insert(std::pair<std::string, XYZType>("IWFM", XYZType::IWFM));
	    it = xyztMap.find(xyzt);
	    if (it != xyztMap.end())
	        return it->second;
	    else {
	        return XYZType::INVALID;
	    }
	}

	/*! \brief The main class that defines all the options
	 *
	 */
	class options {
	public:
	    /*!
	     * The constructor of the option class
	     * @param world_in is the mpi communicator
	     */
		options(boost::mpi::communicator& world_in);

		/*!
		 * Reads the input arguments to the main program.
		 *
		 * See details explanation about the inputs [here](https://docs.microsoft.com/en-us/cpp/cpp/main-function-command-line-args?view=msvc-160)
		 *
		 * @param argc The number of input arguments
		 * @param argv The input arguments themselves
		 * @return true if the reading was successful
		 */
		bool readInput(int argc, char* argv[]);

		/*!
		 *  A getter function to the name of velocity configuration file
		 * @return the velocity configuration file name
		 */
		std::string getVelFname() {
			return velocityFieldFileName;
		}


		ParticleOptions Popt; /// A container of the particle tracking related options
		DomainOptions Dopt; /// A container of the domain related options
		bool gatherMode = false; /// Indicates whether the program should rather particles from files or do particle tracking
		int nproc = 0; /// Number of processors. This should set automatically.
		int niter = 0; /// The number of iterations during the gather mode. This is set also from the input files

		VelType velocityFieldType; /// The velocity field type
		XYZType xyztype;

		std::string Version;
		

		

		
	private:
		boost::mpi::communicator world; /// A reference to the mpi communicator
		// bool bIsMultiThreaded = true;
		//int nThreads = 1;
		std::string velocityFieldFileName;
	};

	options::options(boost::mpi::communicator& world_in)
		:
		world(world_in)
	{
        Version = "0.2.10";
		//if (world.size() > 1) {
		//	bIsMultiThreaded = false;
		//}
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
				std::cout << "| Version : " << Version <<" |" << std::endl;
				std::cout << "|    by  giorgk    |" << std::endl;
				std::cout << "|------------------|" << std::endl;
			}
			return false;
		}

		// Configuration file options
		po::options_description config_options("Configuration file options");
		config_options.add_options()
			// Velocity Options
			("Velocity.XYZType", po::value<std::string>(), "Type of point coordinates e.g CLOUD/IWFM")
			("Velocity.Type", po::value<std::string>(), "Type of velocity. (STEADY or TRANS)")
			("Velocity.ConfigFile", po::value<std::string >(), "Set configuration file for the velocity field")

			// Domain options
			("Domain.Outline", po::value<std::string >(), "A filename that contains the vertices of the outline polygon")
			("Domain.TopFile", po::value<std::string >(), "A filename with the point cloud of the top elevation")
			("Domain.TopRadius", po::value<double>()->default_value(1000), "Search Radius for top elevation")
			("Domain.TopPower", po::value<double>()->default_value(3), "Search Power for top elevation")
			("Domain.BottomFile", po::value<std::string >(), "A filename with the point cloud of the bottom elevation")
			("Domain.BottomRadius", po::value<double>()->default_value(1000), "Search Radius for bottom elevation")
			("Domain.BottomPower", po::value<double>()->default_value(3), "Search Power for bottom elevation")
			("Domain.ProcessorPolys", po::value<std::string >(), "A filename that contains the coordinates of each processor polygon")
			("Domain.ExpandedPolys", po::value<std::string >(), "A filename that contains the coordinates of each Expanded polygon")

			// Stopping criteria
			("StoppingCriteria.MaxIterationsPerStreamline", po::value<int>()->default_value(1000), "Maximum number of steps per streamline")
			("StoppingCriteria.MaxProcessorExchanges", po::value<int>()->default_value(50), "Maximum number of that a particles are allowed to change processors")
			("StoppingCriteria.AgeLimit", po::value<double>()->default_value(-1), "If the particle exceed this limit stop particle tracking. Negative means no limit")
			("StoppingCriteria.StuckIter", po::value<int>()->default_value(10), "After StuckIter exit particle tracking")
            ("StoppingCriteria.AttractFile", po::value<std::string >(), "A filename that contains the particle attractors")
            ("StoppingCriteria.AttractRadius", po::value<double>()->default_value(50), "Particles closer to this distance will be captured from attractors")

			// Step configuration
			("StepConfig.Method", po::value<std::string >(), "Method for stepping")
			("StepConfig.Direction", po::value<double>()->default_value(1), "Backward or forward particle tracking")
			("StepConfig.StepSize", po::value<double>()->default_value(1), "Step Size in units of length")
            ("StepConfig.StepSizeTime", po::value<double>()->default_value(10000), "Step Size in units of time")
            ("StepConfig.nSteps", po::value<int>()->default_value(4), "The number of steps to take within the BBox or the element")
            ("StepConfig.nStepsTime", po::value<int>()->default_value(2), "The number of steps to take within a time step")
			("StepConfig.minExitStepSize", po::value<double>()->default_value(0.1), "Minimum Step Size at the exit")
            //("StepConfig.UpdateStepSize", po::value<int>()->default_value(1), "Update step size wrt bbox")

			// Adaptive Step configurations

			("AdaptStep.MaxStepSize", po::value<double>()->default_value(2), "Maximum Step Size in units of length")
			("AdaptStep.MinStepSize", po::value<double>()->default_value(0.1), "Minimum Step Size in units of length")
			("AdaptStep.increaseRateChange", po::value<double>()->default_value(1.5), "Maximum Step Size in units of length")
			("AdaptStep.limitUpperDecreaseStep", po::value<double>()->default_value(0.75), "Upper limit of decrease step size")
			("AdaptStep.Tolerance", po::value<double>()->default_value(0.1), "Tolerance when the RK45 is used")

			// InputOutput
			("InputOutput.ParticleFile", po::value<std::string >(), "A filename with the initial positions of particles")
			("InputOutput.WellFile", po::value<std::string >(), "A filename with the well locations")
			("InputOutput.OutputFile", po::value<std::string >(), "Prefix for the output file")
			("InputOutput.ParticlesInParallel", po::value<int>()->default_value(1000), "Maximum number run in parallel")
			("InputOutput.GatherOneFile", po::value<int>()->default_value(1), "Put all streamlines into one file")

			// Other
			("Other.Nrealizations", po::value<int>()->default_value(1), "Number of realizations")
			("Other.Version", po::value<std::string >(), "The version of the Ichnos. (Check ichnos.exe -v)")

			// Unsued
			//("Unused.nThreads", po::value<int>()->default_value(1), "Number of threads")
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
				std::cout << "--> Configuration file: " << vm_cmd["config"].as<std::string>().c_str() << std::endl;
				po::store(po::parse_config_file<char>(vm_cmd["config"].as<std::string>().c_str(), config_options), vm_cfg);

                {// Other
                    if (!vm_cfg.count("Other.Version")){
                        std::cout << "Version is not specified!. The current version is " << Version << std::endl;
                        return false;
                    }
                    std::string user_version = vm_cfg["Other.Version"].as<std::string>();
                    if (Version.compare(user_version) != 0){
                        std::cout << "Mismatch between code version (" << Version << ") and Input file version (" << user_version << ")" << std::endl;
                        return false;
                    }
                    Popt.Nrealizations = vm_cfg["Other.Nrealizations"].as<int>();
                }

				{// Velocity Options
					velocityFieldType = castVelType2Enum(vm_cfg["Velocity.Type"].as<std::string>());
					if (velocityFieldType == VelType::INVALID) {
						std::cout << vm_cfg["Velocity.Type"].as<std::string>() << " is an invalid velocity type" << std::endl;
						return false;
					}
					else if (velocityFieldType == VelType::STEADY){// Steady flow field is essentially an one-step transient. The only difference is the input file for XYZ and velocity format
                        velocityFieldType = VelType::TRANS;
					}

					xyztype = castXYZType2Enum(vm_cfg["Velocity.XYZType"].as<std::string>());
					if (xyztype == XYZType::INVALID){
					    std::cout << vm_cfg["Velocity.XYZType"].as<std::string>() << " is an invalid XYZ type" << std::endl;
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
				}

				{// Stopping criteria
					Popt.MaxIterationsPerStreamline = vm_cfg["StoppingCriteria.MaxIterationsPerStreamline"].as<int>();
					Popt.MaxProcessorExchanges = vm_cfg["StoppingCriteria.MaxProcessorExchanges"].as<int>();
					Popt.AgeLimit = vm_cfg["StoppingCriteria.AgeLimit"].as<double>();
					Popt.StuckIterations = vm_cfg["StoppingCriteria.StuckIter"].as<int>();
                    Dopt.AttractorsFile = vm_cfg["StoppingCriteria.AttractFile"].as<std::string>();
                    Dopt.AttractRadius = vm_cfg["StoppingCriteria.AttractRadius"].as<double>();
                    Dopt.AttractRadius = Dopt.AttractRadius * Dopt.AttractRadius;
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

					Popt.StepOpt.dir = Popt.Direction;
					Popt.StepOpt.StepSize = vm_cfg["StepConfig.StepSize"].as<double>();
					Popt.StepOpt.StepSizeTime = vm_cfg["StepConfig.StepSizeTime"].as<double>();
					Popt.StepOpt.nSteps = std::max(vm_cfg["StepConfig.nSteps"].as<int>(),1);
                    Popt.StepOpt.nStepsTime = vm_cfg["StepConfig.nStepsTime"].as<int>();
					Popt.StepOpt.minExitStepSize = vm_cfg["StepConfig.minExitStepSize"].as<double>();
					if (Popt.StepOpt.minExitStepSize < 0 || Popt.StepOpt.minExitStepSize > 1) {
						std::cout << "minExitStepSize should be between 0 and 1. It gets the default value of 0.1" << std::endl;
						Popt.StepOpt.minExitStepSize = 0.1;
					}
                    //Popt.UpdateStepSize = vm_cfg["StepConfig.UpdateStepSize"].as<int>();
                    Popt.UpdateStepSize = 1;
					if (Popt.method == SolutionMethods::RK45){
                        Popt.UpdateStepSize = 0;
					}
				}

				{// Adaptive Step configurations

					Popt.AdaptOpt.MaxStepSize = vm_cfg["AdaptStep.MaxStepSize"].as<double>();
					Popt.AdaptOpt.MinStepSize = vm_cfg["AdaptStep.MinStepSize"].as<double>();
					Popt.AdaptOpt.increaseRateChange = vm_cfg["AdaptStep.increaseRateChange"].as<double>();
					if (Popt.AdaptOpt.increaseRateChange < 1) {
						std::cout << "increaseRateChange should be higher than 1. It gets the default value of 1.5" << std::endl;
						Popt.AdaptOpt.increaseRateChange = 1.5;
					}
					Popt.AdaptOpt.limitUpperDecreaseStep = vm_cfg["AdaptStep.limitUpperDecreaseStep"].as<double>();
					if (Popt.AdaptOpt.limitUpperDecreaseStep < 0 || Popt.AdaptOpt.limitUpperDecreaseStep > 1) {
						std::cout << "limitUpperDecreaseStep should be between 0 and 1. It gets the default value of 0.75" << std::endl;
						Popt.AdaptOpt.limitUpperDecreaseStep = 0.75;
					}
					Popt.AdaptOpt.ToleranceStepSize = vm_cfg["AdaptStep.Tolerance"].as<double>();
				}

				{// InputOutput
					Popt.ParticleFile = vm_cfg["InputOutput.ParticleFile"].as<std::string>();
					Popt.WellFile = vm_cfg["InputOutput.WellFile"].as<std::string>();
					Popt.OutputFile = vm_cfg["InputOutput.OutputFile"].as<std::string>();
					Popt.ParticlesInParallel = vm_cfg["InputOutput.ParticlesInParallel"].as<int>(); 
				}
				//{// Unsued
					//nThreads = vm_cfg["Unused.nThreads"].as<int>();
				//}
			
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
