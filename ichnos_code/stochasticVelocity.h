#pragma once
#include <boost/mpi.hpp>
#include <boost/program_options.hpp>

#include "ichnos_structures.h"
#include "velocity_base.h"

namespace po = boost::program_options;
namespace ic = ICHNOS;

namespace STOCH {

	class MarkovChainVel : public ICHNOS::velocityField {
	public:
		MarkovChainVel(boost::mpi::communicator& world_in);
		void readVelocityField(std::string vf_file, int nPnts);
		void calcVelocity(ic::vec3& vel,
                          std::vector<int>& ids,
                          std::vector<double>& weights,
                          double time = 0);
		void reset();
		void updateStep(double& step);
	private:
		ic::search_tree_stoch Tree;
		// Number of finite states
		int Nstates = -9;
		int Nperiods = -9;
		// A list of transition probability matrices
		ic::TransitionProbabilityMatrix TPM;
		ICHNOS::SingletonGenerator* RG = RG->getInstance();

		bool readNodesVelocities(std::string prefix, std::string suffix,
			std::vector<ic::cgal_point_3>& pp,
			std::vector<ic::STOCH_data>& dd, int leadZeros);

		ic::interpType porType;
		double porosityValue = 1.0;

		bool bIsInitialized = false;
		int currentState;
		int currentPeriod;
		double initial_diameter = 640;
		double initial_ratio = 20;
		int initial_state = 0;
		int initial_period = 0;
		int updateYearType = 0;
		double Scale = 1.0;
		double Threshold;
		double diameter;
		double ratio;
		double Power;
		double search_mult = 2.5;
		double TimePerPeriod = 30;
		double TimeSinceLastUpdate = 0;
		double StepSize;
		ic::vec3 ll, uu, pp, vv;

	};

	MarkovChainVel::MarkovChainVel(boost::mpi::communicator& world_in)
		:
		velocityField(world_in)
	{}

	void MarkovChainVel::reset() {
		bIsInitialized = false;
		diameter = initial_diameter;
		ratio = initial_ratio;
		currentState = initial_state;
		currentPeriod = initial_period;
		TimeSinceLastUpdate = 0;
	}

	void MarkovChainVel::readVelocityField(std::string vf_file, int nPnts) {
        if (world.rank() == 0)
		    std::cout << "Velocity configuration file: " << vf_file << std::endl;
		po::options_description velocityFieldOptions("Velocity field options");
		po::variables_map vm_vfo;
		velocityFieldOptions.add_options()
			("Prefix", po::value<std::string>(), "Prefix for the input files.")
			("Suffix", po::value<std::string>(), "Suffix for the input files.")
			("OwnerThreshold", po::value<double>()->default_value(0.75), "Threshold for the processor ownership")
			("VelocityMultiplier", po::value<double>()->default_value(1), "This is a multiplier to scale velocity")
			("LeadingZeros", po::value<int>()->default_value(4), "e.g 0002->4, 000->3")
			("Nstates", po::value<int>(), "Number of states of the Markov chain model")
			("NPeriods", po::value<int>(), "Number of periods of the Markov chain model")
			("TPmatrix", po::value<std::string>(), "Transition probability matrix")
			("Porosity", po::value<std::string>(), "Porocity. Either a file or a single number")
			("InitState", po::value<int>()->default_value(1), "An integer that represents the initial state [0-Nstates)")//The input must be based on 1. It is converted to zero based in code
			("InitPeriod", po::value<int>()->default_value(1), "An integer that represents the initial period [0-NPeriods)")//The input must be based on 1. It is converted to zero based in code
			("UpdateYearType", po::value<int>()->default_value(1), "An integer that represents Which period the water year type should be updated")//The input must be based on 1. It is converted to zero based in code
			("StepSize", po::value<double>()->default_value(1), "Step size in time units")// If the velocity is m/day then the step is in days
			("TimePerPeriod", po::value<double>()->default_value(30), "The time of a period.")// if the step size is days and periods months this should be 30.
			("InitDiameter", po::value<double>()->default_value(5000), "Initial diameter")
			("InitRatio", po::value<double>()->default_value(1), "Initial ratio")
			("Power", po::value<double>()->default_value(3.0), "Power of the IDW interpolation")
			("Scale", po::value<double>()->default_value(1.0), "Scale the domain before velocity calculation")
			;

		po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions), vm_vfo);

		// Transition probabilities
		Nstates = vm_vfo["Nstates"].as<int>();
		initial_state = vm_vfo["InitState"].as<int>() - 1;
		if (initial_state >= Nstates)
			initial_state = 0;
		Nperiods = vm_vfo["NPeriods"].as<int>();
		initial_period = vm_vfo["InitPeriod"].as<int>() - 1;
		updateYearType = vm_vfo["UpdateYearType"].as<int>() - 1;
		if (initial_period >= Nperiods)
			initial_period = 0;
		
		if (Nstates > 0) {
			TPM.setNstates(Nstates);
			TPM.readData(vm_vfo["TPmatrix"].as<std::string>());
		}
		else {
			std::cerr << "The number of states must be positive" << std::endl;
		}

		Power = vm_vfo["Power"].as<double>();
		Scale = vm_vfo["Scale"].as<double>();
		Nperiods = vm_vfo["NPeriods"].as<int>();
		TimePerPeriod = vm_vfo["TimePerPeriod"].as<double>();
		initial_diameter = vm_vfo["InitDiameter"].as<double>();
		initial_ratio = vm_vfo["InitRatio"].as<double>();
		StepSize = vm_vfo["StepSize"].as<double>();
		//std::string TPmatrixFile = vm_vfo["TPmatrix"].as<std::string>();

		//bool tf = readTransitionProbabilities(TPmatrixFile);

		std::string prefix = vm_vfo["Prefix"].as<std::string>();
		std::string suffix = vm_vfo["Suffix"].as<std::string>();
		int leadZeros = vm_vfo["LeadingZeros"].as<int>();

		{ // Read Velocities
			std::vector<ic::cgal_point_3> pp;
			std::vector<ic::STOCH_data> dd;
			bool tf = readNodesVelocities(prefix, suffix, pp, dd, leadZeros);

			{// Build tree
				auto start = std::chrono::high_resolution_clock::now();
				Tree.insert(boost::make_zip_iterator(boost::make_tuple(pp.begin(), dd.begin())),
					boost::make_zip_iterator(boost::make_tuple(pp.end(), dd.end())));
				Tree.build();
				auto finish = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = finish - start;
				std::cout << "Point Set Building time: " << elapsed.count() << std::endl;
			}
		}

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
					// This is not implemented yet
					std::cout << "Variable porocity not implemented yet. Setting porocity equal to 1" << std::endl;
					porType = ic::interpType::SCALAR;
                    porosityValue = 1.0;
				}
			}
		}
		else {
			porType = ic::interpType::INGORE;
            porosityValue = 1.0;
		}

	}

	void MarkovChainVel::calcVelocity(ic::vec3& vel,
                                      std::vector<int>& ids,
                                      std::vector<double>& weights,
                                      double time) {
	    ic::vec3 p;
		// If this is the first point of this streamline we will carry out one additional range search
		ll.zero();
		uu.zero();
		pp.zero();
		vv.zero();
		if (!bIsInitialized) {
			ic::calculate_search_box(p, ll, uu, diameter, ratio, search_mult);
			ic::cgal_point_3 llp(ll.x, ll.y, ll.z);
			ic::cgal_point_3 uup(uu.x, uu.y, uu.z);
			ic::Fuzzy_iso_box_stoch fib(llp, uup, 0.0);
			std::vector<boost::tuples::tuple<ic::cgal_point_3, ic::STOCH_data>> tmp;
			Tree.search(std::back_inserter(tmp), fib);
			if (tmp.size() == 0) {
				vel = ic::vec3(-99999, -99999, -99999);
				return;
			}
			// Find the closest point
			double mindist = 1000000000;
			double tmp_diam;
			double tmp_ratio;
			ic::STOCH_data closest_point_data;
			for (unsigned int i = 0; i < tmp.size(); ++i) {
				double dist = p.distance(tmp[i].get<0>().x(), tmp[i].get<0>().y(), tmp[i].get<0>().z());
				if (dist < mindist) {
					mindist = dist;
					tmp_diam = tmp[i].get<1>().diameter;
					tmp_ratio = tmp[i].get<1>().ratio;
				}
			}
			diameter = tmp_diam;
			ratio = tmp_ratio;
			bIsInitialized = true;
		}

		{
			auto start = std::chrono::high_resolution_clock::now();
			std::map<int, double>::iterator itd;
			std::vector<boost::tuples::tuple<ic::cgal_point_3, ic::STOCH_data>> tmp;
			while (true) {
				tmp.clear();
				ic::calculate_search_box(p, ll, uu, diameter, ratio, search_mult);
				ic::cgal_point_3 llp(ll.x, ll.y, ll.z);
				ic::cgal_point_3 uup(uu.x, uu.y, uu.z);
				ic::Fuzzy_iso_box_stoch fib(llp, uup, 0.0);
				Tree.search(std::back_inserter(tmp), fib);
				if (tmp.size() >= 3) {
					break;
				}
				else {
					diameter = diameter * 1.5;
					if (diameter > initial_diameter) {
						vel = ic::vec3(-99999, -99999, -99999);
						return;
					}
				}
			}

			// por_xyz refers to point translated to origin
			double porx, pory, porz, scaled_dist, actual_dist, w;
			bool calc_average = true;
			double sumW = 0;
			ic::vec3 sumWVal;
			ic::vec3 val;
			double mindist = 1000000000;
			double tmp_diam;
			double tmp_ratio;

			for (unsigned int i = 0; i < tmp.size(); ++i) {
				porx = p.x - tmp[i].get<0>().x();
				pory = p.y - tmp[i].get<0>().y();
				porz = p.z - tmp[i].get<0>().z();
				actual_dist = std::sqrt(porx * porx + pory * pory + porz * porz);
				if (actual_dist < mindist) {
					mindist = actual_dist;
					tmp_diam = tmp[i].get<1>().diameter;
					tmp_ratio = tmp[i].get<1>().ratio;
				}
				porz = porz * ratio * Scale + porz * (1 - Scale);
				scaled_dist = std::sqrt(porx * porx + pory * pory + porz * porz);

				if (actual_dist < Threshold) {
					int n = tmp[i].get<1>().v.nValues(currentState, currentPeriod);
					int rindex = RG->randomNumber(0, n);
					tmp[i].get<1>().v.getValue(currentState, currentPeriod, rindex, vel);
					calc_average = false;
				}
				else {
					w = 1 / std::pow(scaled_dist, Power);
//					itd = proc_map.find(tmp[i].get<1>().proc);
//					if (itd == proc_map.end()) {
//						proc_map.insert(std::pair<int, double>(tmp[i].get<1>().proc, w));
//					}
//					else {
//						itd->second += w;
//					}
					sumW += w;
					int n = tmp[i].get<1>().v.nValues(currentState, currentPeriod);
					int rindex = RG->randomNumber(0, n);
					tmp[i].get<1>().v.getValue(currentState, currentPeriod, rindex, val);
					sumWVal = sumWVal + val * w;
				}
			}
			if (tmp.size() > 50) { // Double check that and explain why
				diameter = tmp_diam;
				ratio = tmp_ratio;
			}

//			itd = proc_map.begin();
//			for (; itd != proc_map.end(); ++itd) {
//				itd->second = itd->second / sumW;
//			}
			if (calc_average)
				vel = sumWVal * (1 / sumW);
			


			auto finish = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = finish - start;
		}

		if (porType == ic::interpType::CLOUD) {
			// TODO not implemented yet
		}
		else if (porType == ic::interpType::SCALAR) {
			vel = vel * (1 / porosityValue);
		}

        // TODO Understand why this is needed
		pp = p; // I dont understand this. It seems we need to use these variables in the updateStep method
		vv = vel;
	}

	struct tempND {
		int proc;
		double x;
		double y;
		double diameter;
		std::vector<double> elev;
		std::vector<double> ratio;
		std::vector<ic::Stochastic_Velocity> vel;
	};

	bool MarkovChainVel::readNodesVelocities(std::string prefix, std::string suffix,
		std::vector<ic::cgal_point_3>& pp,
		std::vector<ic::STOCH_data>& dd, int leadZeros) {

		auto start = std::chrono::high_resolution_clock::now();
		std::string filexy = prefix + "XYZ_" + ic::num2Padstr(/*dbg_rank*/world.rank(), leadZeros) + suffix;
		std::string filevel = prefix + "VEL_" + ic::num2Padstr(/*dbg_rank*/world.rank(), leadZeros) + suffix;

		std::map<int, tempND> ND;
		std::map<int, tempND>::iterator itND;


		// -------------READ XY Velocity NODE information --------------------------
		std::cout << "Proc " << world.rank() << " is reading the XYZ nodes" << std::endl;
		std::ifstream datafileXY(filexy.c_str());
		if (!datafileXY.good()) {
			std::cout << "Can't open the file " << filexy << std::endl;
		}
		else {
			std::string line;
			int id, nlay;
			double u;
			ic::Stochastic_Velocity empty_stoch_vel(Nstates, Nperiods);
			empty_stoch_vel.initStatesPeriods();

			

			while (getline(datafileXY, line)) {
				if (line.size() > 1) {
					std::istringstream inp(line.c_str());
					tempND tmp;
					inp >> id;
					inp >> tmp.proc;
					inp >> tmp.diameter;
					inp >> tmp.x;
					inp >> tmp.y;
					inp >> nlay;
					tmp.elev.resize(nlay);
					tmp.ratio.resize(nlay);
					//tmp.vel.resize(nlay);
					for (int i = 0; i < nlay; i++) {
						inp >> u;
						tmp.elev[i] = u;
						tmp.vel.push_back(empty_stoch_vel);
					}
					/*
					for (int i = 0; i < nlay; i++) {
						if (i == 0) {
							double dz = tmp.elev[0] - tmp.elev[1];
							tmp.ratio[i] = tmp.diameter / dz;
						}
						else if (i == nlay - 1) {
							double dz = tmp.elev[i-1] - tmp.elev[i];
							tmp.ratio[i] = tmp.diameter / dz;
						}
						else {
							double dz = (tmp.elev[i - 1] + tmp.elev[i]) / 2 - (tmp.elev[i] + tmp.elev[i + 1]) / 2;
							tmp.ratio[i] = tmp.diameter / dz;
						}
					}
					*/
					ND.insert(std::pair<int, tempND>(id, tmp));
				}
			}
			datafileXY.close();
		}

		// -------------READ Velocity information of the nodes -------------------
		
		std::cout << "Proc " << world.rank() << " is reading the Velocity of the nodes" << std::endl;
		std::ifstream datafileVel(filevel.c_str());
		if (!datafileVel.good()) {
			std::cout << "Can't open the file " << filevel << std::endl;
		}
		else {
			std::string line;
			int id, lay, istate, imon, nVel;
			double u, ratio;
			while (getline(datafileVel, line)) {

				if (line.size() > 1) {
					std::istringstream inp(line.c_str());
					inp >> id;
					inp >> lay;
					inp >> istate;
					inp >> imon;
					inp >> ratio;
					inp >> nVel;
					std::vector<ic::vec3> vel_vec(nVel, ic::vec3());
					itND = ND.find(id);
					itND->second.ratio[lay - 1] = ratio;
					for (int j = 0; j < 3; j++) {
						for (int i = 0; i < nVel; i++) {
							inp >> u;
							if (j == 0)
								vel_vec[i].x = u;
							else if (j == 1)
								vel_vec[i].y = u;
							else if (j == 2)
								vel_vec[i].z = u;
						}
					}
					for (int i = 0; i < nVel; i++) {
						itND->second.vel[lay-1].addValue(istate-1, imon-1, vel_vec[i]);
					}
				}
			}
			datafileVel.close();
			std::cout << "Proc " << world.rank() << " finish reading the Velocity of the nodes" << std::endl;
		}

		int cnt_id = 0;
		for (itND = ND.begin(); itND != ND.end(); itND++) {
			for (unsigned int i = 0; i < itND->second.elev.size(); i++) {
				pp.push_back(ic::cgal_point_3(itND->second.x, itND->second.y, itND->second.elev[i]));
				ic::STOCH_data tmp;
				tmp.id = cnt_id;
				cnt_id++;
				tmp.diameter = itND->second.diameter;
				tmp.proc = itND->second.proc;
				tmp.ratio = itND->second.ratio[i];
				tmp.v = itND->second.vel[i];
				dd.push_back(tmp);
			}
		}
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Read data time: " << elapsed.count() << std::endl;
		return true;
	}


	void MarkovChainVel::updateStep(double& step) {
		TimeSinceLastUpdate += StepSize;
		step = vv.len() * StepSize;
		if (TimeSinceLastUpdate >= TimePerPeriod) {
			TimeSinceLastUpdate = 0;
			currentPeriod++;
			if (currentPeriod >= Nperiods)
				currentPeriod = 0;
		}
		if (currentPeriod == updateYearType) {
			double r = RG->randomNumber();
			currentState = TPM.nextState(currentState, pp, r);
		}
	}
}