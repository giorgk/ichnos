#pragma once
#include <boost/mpi.hpp>
#include <boost/program_options.hpp>

#include "ichnos_structures.h"
#include "velocity_base.h"
#include "ichnos_XYZ_base.h"
#include "TransientVelocity.h"
#include "velocity_mesh2D.h"

namespace po = boost::program_options;
namespace ic = ICHNOS;

namespace ICHNOS {

	class MarkovChainVelCloud : public ICHNOS::CloudVel {
	public:
		MarkovChainVelCloud(boost::mpi::communicator& world_in, ic::XYZ_base &XYZ_in);
		bool readVelocityField(std::string vf_file);
		void calcVelocity(ic::vec3& p, ic::vec3& vel,
                          std::map<int, double>& proc_map,
                          ic::helpVars& pvlu,
                          bool& out,
                          double time = 0);
		void reset(Streamline& S);
		void updateStep(helpVars& pvlu);
	private:
		//ic::search_tree_stoch Tree;
		// Number of finite states
		int Nstates = -9;
		int Nperiods = -9;
		// A list of transition probability matrices
		ic::TransitionProbabilityMatrix TPM;
		ICHNOS::SingletonGenerator* RG = RG->getInstance();
        MCVpools MCVP;

		//bool readNodesVelocities(std::string prefix, std::string suffix,
		//	std::vector<ic::cgal_point_3>& pp,
		//	std::vector<ic::STOCH_data>& dd, int leadZeros);

		bool bIsInitialized = false;
		int currentState;
		int currentPeriod;
		int initial_state = 0;
		int initial_period = 0;
		int UpdatePeriod = 0;
		//double Scale = 1.0;
		//double search_mult = 2.5;
		//double TimePerPeriod = 30;
		double TimeSinceLastUpdate = 0;
		//double StepSize;
		//ic::vec3 ll, uu, pp, vv;
		void updatePeriod(bool add);

	};

	MarkovChainVelCloud::MarkovChainVelCloud(boost::mpi::communicator& world_in, ic::XYZ_base &XYZ_in)
		:
        CloudVel(world_in, XYZ_in)
	{}

	void MarkovChainVelCloud::reset(Streamline& S) {
		CloudVel::reset(S);
		bIsInitialized = false;
		//diameter = initial_diameter;
		//ratio = initial_ratio;
		currentState = initial_state;
		currentPeriod = initial_period;
		TimeSinceLastUpdate = 0;
	}

	bool MarkovChainVelCloud::readVelocityField(std::string vf_file) {
        CloudVel::readVelocityField(vf_file);
        if (world.rank() == 0)
		    std::cout << "Markov Chain related parameters: " << vf_file << std::endl;

		po::options_description velocityFieldOptions("Velocity field options");
		po::variables_map vm_vfo;
		velocityFieldOptions.add_options()
			// MarkovChain
			("MarkovChain.Nstates", po::value<int>(), "Number of states of the Markov chain model")
			("MarkovChain.Nperiods", po::value<int>(), "Number of periods of the Markov chain model")
			("MarkovChain.UpdatePeriod", po::value<int>()->default_value(1), "An integer that shows after which period the state should be updated")//The input must be based on 1. It is converted to zero based in code
            ("MarkovChain.InitState", po::value<int>()->default_value(1), "The initial State for each particle")
            ("MarkovChain.InitPeriod", po::value<int>()->default_value(1), "The initial Period for each particle")
            ("MarkovChain.TPmatrix", po::value<std::string>(), "Transition probability matrix")
            ("MarkovChain.StateSequence", po::value<std::string>(), "State Sequence file")
            ("MarkovChain.VelocityPools", po::value<std::string>(), "Velocity pools")
			("MarkovChain.TimesPerPeriod", po::value<std::string>(), "Time per period. it is used with TPmatrix option only")
			;

		po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions, true), vm_vfo);

		// Transition probabilities
		Nstates = vm_vfo["MarkovChain.Nstates"].as<int>();
		initial_state = vm_vfo["MarkovChain.InitState"].as<int>() - 1;
		if (initial_state >= Nstates)
			initial_state = 0;
		Nperiods = vm_vfo["MarkovChain.Nperiods"].as<int>();
		initial_period = vm_vfo["MarkovChain.InitPeriod"].as<int>() - 1;
		if (initial_period >= Nperiods)
			initial_period = 0;
		UpdatePeriod = vm_vfo["MarkovChain.UpdatePeriod"].as<int>() - 1;

		
		if (Nstates > 0) {
			TPM.setNstates(Nstates);
            std::string tpmatrixfile = vm_vfo["MarkovChain.TPmatrix"].as<std::string>();
            if (!tpmatrixfile.empty()){
                TPM.readData(tpmatrixfile);

				std::string timesperiodfile = vm_vfo["MarkovChain.TimesPerPeriod"].as<std::string>();
				TPM.readTimesPerPeriod(timesperiodfile);
            }
            else{
                std::string sequenfile = vm_vfo["MarkovChain.StateSequence"].as<std::string>();
                if (sequenfile.empty()){
                    std::cout << "Both TPmatrix and StateSequence are empty!!" << std::endl;
                    return false;
                }
                else{
                    //TODO
                    //Write reader for reading sequence file
                }
            }

            MCVP.init(Nstates, Nperiods);
            std::string poolfile = vm_vfo["MarkovChain.VelocityPools"].as<std::string>();
            bool tf = MCVP.readPoolData(poolfile);
            return tf;
		}
		else {
			std::cout << "The number of states must be positive" << std::endl;
            return false;
		}
        return false;
	}

	void MarkovChainVelCloud::calcVelocity(ic::vec3& p, ic::vec3& vel,
                                           std::map<int, double>& proc_map,
                                           ic::helpVars& pvlu,
                                           bool& out,
                                           double tm) {
        std::vector<int> ids;
        std::vector<double> weights;
        XYZ.calcWeights(p, ids,weights, proc_map, pvlu, out);
        if (!out)
            return;

        // Pick a random time index from the pool
        int poolSize = MCVP.poolSize(currentState, currentPeriod);
        int idx = RG->randomNumber(0, poolSize);
        int itime = MCVP.getVelocity(currentState, currentPeriod,idx);

        double sumW = 0;
        ic::vec3 sumWVal;
        vel.zero();
        for (unsigned int i = 0; i < ids.size(); ++i) {
            vel = VEL.getVelocity(ids[i], itime, itime, 0.0);
            sumW += weights[i];
            sumWVal = sumWVal + vel * weights[i];
        }
        vel = sumWVal * (1 / sumW);

        double por = porosity.interpolate(p);
        vel = vel * (1/por);

        /*
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
        */
	}

	//struct tempND {
	//	int proc;
	//	double x;
	//	double y;
	//	double diameter;
	//	std::vector<double> elev;
	//	std::vector<double> ratio;
	//	std::vector<ic::Stochastic_Velocity> vel;
	//};

    /*
	bool MarkovChainVelCloud::readNodesVelocities(std::string prefix, std::string suffix,
                                                  std::vector<ic::cgal_point_3>& pp,
                                                  std::vector<ic::STOCH_data>& dd, int leadZeros) {

		auto start = std::chrono::high_resolution_clock::now();
		std::string filexy = prefix + "XYZ_" + ic::num2Padstr(world.rank(), leadZeros) + suffix;
		std::string filevel = prefix + "VEL_" + ic::num2Padstr(world.rank(), leadZeros) + suffix;

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
					/
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
					/
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
    */


	void MarkovChainVelCloud::updateStep(helpVars& pvlu) {
		double currentPeriodLength = TPM.periodLength(currentPeriod);
		if (stepOpt.dir > 0){
			TimeSinceLastUpdate +=pvlu.actualStep;
			while (true){
				if (TimeSinceLastUpdate < currentPeriodLength){
					break;
				}
				TimeSinceLastUpdate = TimeSinceLastUpdate - currentPeriodLength;
				updatePeriod(true);
				currentPeriodLength = TPM.periodLength(currentPeriod);
				if (currentPeriod == UpdatePeriod) {
					double r = RG->randomNumber();
					currentState = TPM.nextState(currentState, pvlu.pp, r);
				}

			}
		}
		else{
			TimeSinceLastUpdate -=pvlu.actualStep;
			while (true){
				if (std::abs(TimeSinceLastUpdate) < currentPeriodLength){
					break;
				}
				TimeSinceLastUpdate = TimeSinceLastUpdate + currentPeriodLength;
				updatePeriod(false);
				currentPeriodLength = TPM.periodLength(currentPeriod);
				if (currentPeriod == UpdatePeriod) {
					double r = RG->randomNumber();
					currentState = TPM.nextState(currentState, pvlu.pp, r);
				}
			}
		}
	}

	void MarkovChainVelCloud::updatePeriod(bool add) {
		if (add){
			currentPeriod += 1;
		}
		else{
			currentPeriod -= 1;
		}
		if (currentPeriod >= Nperiods){
			currentPeriod = 0;
		}
		else if (currentPeriod < 0){
			currentPeriod = Nperiods - 1;
		}
	}



    class MarkovChainVelMESH2D : public ICHNOS::Mesh2DVel{
    public:
        MarkovChainVelMESH2D(boost::mpi::communicator& world_in, ic::XYZ_base &XYZ_in);
        bool readVelocityField(std::string vf_file);
        void calcVelocity(ic::vec3& p, ic::vec3& vel,
                          std::map<int, double>& proc_map,
                          ic::helpVars& pvlu,
                          bool& out,
                          double time = 0);
        void reset(Streamline& S);
        void updateStep(helpVars& pvlu);
    private:
        int Nstates = -9;
        int Nperiods = -9;
        ic::TransitionProbabilityMatrix TPM;
        ICHNOS::SingletonGenerator* RG = RG->getInstance();
        MCVpools MCVP;

        bool bIsInitialized = false;
        int currentState;
        int currentPeriod;
        int initial_state = 0;
        int initial_period = 0;
        int UpdatePeriod = 0;

        double TimeSinceLastUpdate = 0;
        void updatePeriod(bool add);
    };

    MarkovChainVelMESH2D::MarkovChainVelMESH2D(boost::mpi::communicator &world_in, ic::XYZ_base &XYZ_in)
        :
        Mesh2DVel(world_in, XYZ_in)
    {}

    void MarkovChainVelMESH2D::reset(ICHNOS::Streamline &S) {
        Mesh2DVel::reset(S);
        bIsInitialized = false;
        currentState = initial_state;
        currentPeriod = initial_period;
        TimeSinceLastUpdate = 0;
    }

    bool MarkovChainVelMESH2D::readVelocityField(std::string vf_file) {
        Mesh2DVel::readVelocityField(vf_file);
        if (world.rank() == 0)
            std::cout << "Markov Chain related parameters: " << vf_file << std::endl;

        po::options_description velocityFieldOptions("Velocity field options");
        po::variables_map vm_vfo;
        velocityFieldOptions.add_options()
            // MarkovChain
            ("MarkovChain.Nstates", po::value<int>(), "Number of states of the Markov chain model")
            ("MarkovChain.Nperiods", po::value<int>(), "Number of periods of the Markov chain model")
            ("MarkovChain.UpdatePeriod", po::value<int>()->default_value(1), "An integer that shows after which period the state should be updated")//The input must be based on 1. It is converted to zero based in code
            ("MarkovChain.InitState", po::value<int>()->default_value(1), "The initial State for each particle")
            ("MarkovChain.InitPeriod", po::value<int>()->default_value(1), "The initial Period for each particle")
            ("MarkovChain.TPmatrix", po::value<std::string>(), "Transition probability matrix")
            ("MarkovChain.StateSequence", po::value<std::string>(), "State Sequence file")
            ("MarkovChain.VelocityPools", po::value<std::string>(), "Velocity pools")
            ("MarkovChain.TimesPerPeriod", po::value<std::string>(), "Time per period. it is used with TPmatrix option only")
            ;

        po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions, true), vm_vfo);
        // Transition probabilities
        Nstates = vm_vfo["MarkovChain.Nstates"].as<int>();
        initial_state = vm_vfo["MarkovChain.InitState"].as<int>() - 1;
        if (initial_state >= Nstates)
            initial_state = 0;
        Nperiods = vm_vfo["MarkovChain.Nperiods"].as<int>();
        initial_period = vm_vfo["MarkovChain.InitPeriod"].as<int>() - 1;
        if (initial_period >= Nperiods)
            initial_period = 0;
        UpdatePeriod = vm_vfo["MarkovChain.UpdatePeriod"].as<int>() - 1;

        if (Nstates > 0){
            TPM.setNstates(Nstates);
            std::string tpmatrixfile = vm_vfo["MarkovChain.TPmatrix"].as<std::string>();
            if (!tpmatrixfile.empty()){
                TPM.readData(tpmatrixfile);
                std::string timesperiodfile = vm_vfo["MarkovChain.TimesPerPeriod"].as<std::string>();
                TPM.readTimesPerPeriod(timesperiodfile);
            }
            else{
                std::string sequenfile = vm_vfo["MarkovChain.StateSequence"].as<std::string>();
                if (sequenfile.empty()){
                    std::cout << "Both TPmatrix and StateSequence are empty!!" << std::endl;
                    return false;
                }
                else{
                    //TODO
                    //Write reader for reading sequence file
                }
            }

            MCVP.init(Nstates, Nperiods);
            std::string poolfile = vm_vfo["MarkovChain.VelocityPools"].as<std::string>();
            bool tf = MCVP.readPoolData(poolfile);
            return tf;
        }
        else{
            std::cout << "The number of states must be positive" << std::endl;
            return false;
        }
        return false;
    }

    void MarkovChainVelMESH2D::calcVelocity(ic::vec3 &p, ic::vec3 &vel,
                                            std::map<int, double> &proc_map,
                                            ic::helpVars &pvlu,
                                            bool &out,
                                            double time) {

        std::vector<int> ids;
        std::vector<double> weights;
        XYZ.calcWeights(p, ids,weights, proc_map, pvlu, out);
        if (!out){
            vel.makeInvalid();
            return;
        }

        if (ids[0] == -9 || ids[1] == -9){
            vel.makeInvalid();
            return;
        }

        // Pick a random time index from the pool
        int poolSize = MCVP.poolSize(currentState, currentPeriod);
        int idx = RG->randomNumber(0, poolSize);
        int itime = MCVP.getVelocity(currentState, currentPeriod,idx);

        switch (interp_type){
            case ic::MeshVelInterpType::ELEMENT:
            {
                elementInterpolation(vel, ids, itime, itime, 0.0);
                break;
            }
            case ic::MeshVelInterpType::NODE:
            {
                nodeInterpolation(vel, ids, weights, itime, itime, 0.0);
                break;
            }
            case ic::MeshVelInterpType::FACE:
            {
                faceInterpolation(vel, ids, weights, itime, itime, 0.0);
                break;
            }
        }
        double por = porosity.interpolate(p);
        vel = vel * (1/por);
    }

    void MarkovChainVelMESH2D::updateStep(ICHNOS::helpVars &pvlu) {
        double currentPeriodLength = TPM.periodLength(currentPeriod);
        if (stepOpt.dir > 0){
            TimeSinceLastUpdate +=pvlu.actualStep;
            while (true){
                if (TimeSinceLastUpdate < currentPeriodLength){
                    break;
                }
                TimeSinceLastUpdate = TimeSinceLastUpdate - currentPeriodLength;
                updatePeriod(true);
                currentPeriodLength = TPM.periodLength(currentPeriod);
                if (currentPeriod == UpdatePeriod) {
                    double r = RG->randomNumber();
                    currentState = TPM.nextState(currentState, pvlu.pp, r);
                }
            }
        }
        else{
            TimeSinceLastUpdate -=pvlu.actualStep;
            while (true){
                if (std::abs(TimeSinceLastUpdate) < currentPeriodLength){
                    break;
                }
                TimeSinceLastUpdate = TimeSinceLastUpdate + currentPeriodLength;
                updatePeriod(false);
                currentPeriodLength = TPM.periodLength(currentPeriod);
                if (currentPeriod == UpdatePeriod) {
                    double r = RG->randomNumber();
                    currentState = TPM.nextState(currentState, pvlu.pp, r);
                }
            }
        }
    }

    void MarkovChainVelMESH2D::updatePeriod(bool add) {
        if (add){
            currentPeriod += 1;
        }
        else{
            currentPeriod -= 1;
        }
        if (currentPeriod >= Nperiods){
            currentPeriod = 0;
        }
        else if (currentPeriod < 0){
            currentPeriod = Nperiods - 1;
        }
    }



}