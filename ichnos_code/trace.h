#pragma once
#include <boost/mpi.hpp>

#include "ichnos_structures.h"
#include "velocity_base.h"
#include "ichnos_domain.h"
#include "ichnos_XYZ_base.h"


namespace ICHNOS {

	class ParticleTrace {
	public:
		ParticleTrace(
			boost::mpi::communicator& world_in,
            XYZ_base &XYZ_in,
			velocityField &VF_in,
			DomainBase &dom_in,
			ParticleOptions& popt_in);

		void Trace();
		void Trace(int ireal);

		bool readInputFiles(std::vector<Streamline>& S);

	private:
		const ParticleOptions popt;
		velocityField& VF;
		DomainBase& Domain;
        XYZ_base& XYZ;
		boost::mpi::communicator world;
		bool findNextPoint(const Particle& P, vec3& pnew, vec3& vnew, double& tm, ExitReason& er);
		ExitReason traceInner(Streamline& S);
		void traceOuter(std::vector<Streamline>& S, int iter, int ireal);
		ExitReason CheckNewPoint(vec3& p);
		ExitReason CheckNewPointAndCalcVelocity(vec3& p, vec3& v/*, int &proc*/, double tm = 0);

		bool EulerStep(const Particle& P, vec3& pnew, vec3& vnew, double& tm, ExitReason& er);
		bool RK2Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, ExitReason& er);
		bool RK4Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, ExitReason& er);
		bool RK45Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, ExitReason& er);

		void updateStep(double currentAge);

		const std::vector<std::vector<double> > CF = getRK45coef();

		double adaptStepSize;
		double actualStep;
		double actualStepTime;
		vec3 pp, vv, ll, uu;
	};

	ParticleTrace::ParticleTrace(
		boost::mpi::communicator& world_in,
        XYZ_base &XYZ_in,
		velocityField& VF_in,
		DomainBase& dom_in,
		ParticleOptions &popt_in)
		:
		world(world_in),
        XYZ(XYZ_in),
		VF(VF_in),
		Domain(dom_in),
		popt(popt_in)
	{
		adaptStepSize = popt.StepOpt.StepSize;
	}

	bool ParticleTrace::readInputFiles(std::vector <Streamline>& S) {
		bool tfw = false;
		bool tfp = false;
		if (!popt.WellFile.empty()) {
			tfw = READ::readWellFile(popt.WellFile, S, VF.getVelType());
		}

		if (!popt.ParticleFile.empty()) {
			tfp = READ::readParticleFile(popt.ParticleFile, S, VF.getVelType());
		}

		if (tfw || tfp)
			return true;
		else
			return false;
	}

	void ParticleTrace::Trace() {
		Trace(0);
	}

	void ParticleTrace::Trace(int ireal) {
		int my_rank = world.rank();
		int nproc = world.size();
		// Make sure all processors are here
		world.barrier();

		// The first processor reads the well and particle files
		std::vector<Streamline> ALL_Streamlines;
		std::vector<Streamline> Part_Streamlines;
		if (my_rank == 0) {
			bool tf = readInputFiles(ALL_Streamlines);
		}
		// All processors have to wait for the master
		world.barrier();
		ICHNOS::SingletonGenerator* RG = RG->getInstance();

		int outer_iter = 0;
		int particle_index_start = 0;
        int particle_index_end = 0;
		while (true) {
			Part_Streamlines.clear();
			world.barrier();
			if (my_rank == 0) {
			    particle_index_end = particle_index_start + popt.ParticlesInParallel;
			    if (particle_index_end + static_cast<int>(static_cast<double>(popt.ParticlesInParallel)*0.3) > static_cast<int>(ALL_Streamlines.size())){
                    particle_index_end = ALL_Streamlines.size();
			    }
			    for (int i = particle_index_start; i < particle_index_end; ++i){
                    Part_Streamlines.push_back(ALL_Streamlines[i]);
			    }


				/*
				int N = popt.ParticlesInParallel;
				int minPartinparallel = 1000;
				//std::cout << "@#&%$&# Dont forget to harcode the minPartinparallel back to 1000. Currently is " << minPartinparallel << std::endl;
				if (static_cast<int>(ALL_Streamlines.size()) < popt.ParticlesInParallel + minPartinparallel) {
					N = ALL_Streamlines.size();
				}
				for (int i = 0; i < N; ++i) {
					if (ALL_Streamlines.size() == 0)
						break;
					int ii = RG->randomNumber(0, ALL_Streamlines.size());
					Part_Streamlines.push_back(ALL_Streamlines[ii]);
					ALL_Streamlines.erase(ALL_Streamlines.begin() + ii);
				}
				 std::cout << "||	Simulating " << Part_Streamlines.size() << " Streamlines. # streamlines to trace: "  << ALL_Streamlines.size() << std::endl;
				*/
				std::cout << "||    Simulating from " << particle_index_start << " - " << particle_index_end-1 <<" out of " << ALL_Streamlines.size() << std::endl;
			}
			world.barrier();
			//std::cout << my_rank << " Before: " << Part_Streamlines.size() << std::endl;
			MPI::Sent_receive_Initial_streamlines_from0(Part_Streamlines, my_rank, world);
			//std::cout << my_rank << " After: " << Part_Streamlines.size() << std::endl;

			traceOuter(Part_Streamlines, outer_iter, ireal);

            int nRemainingstreamlines = 0;
			if (my_rank == 0){
                nRemainingstreamlines = ALL_Streamlines.size() - particle_index_end;
                particle_index_start = particle_index_end;
			}

			//int nRemainingstreamlines = ALL_Streamlines.size();
			//std::cout << "Proc " << my_rank << " N_part2send = " << N_part2send << std::endl;
			MPI::sumScalar<int>(nRemainingstreamlines, nproc, world, MPI_INT);
			outer_iter++;
			
			if (nRemainingstreamlines == 0)
				break;
			
		}
		if (my_rank == 0 && world.size() > 1 && ireal == popt.Nrealizations - 1) {
			std::cout << "Particle tracking completed. You can now gather the streamlines by running:" << std::endl;
			std::cout << "./ichnos -c " << popt.configfile << " -g -n " << nproc << " -i " << outer_iter << std::endl;
		}
	}

	void ParticleTrace::traceOuter(std::vector<Streamline>& S, int iter, int ireal) {
		int my_rank = /*dbg_rank;*/ world.rank();
		int n_proc = world.size();

		boostPolygon thisDomainPolygon = Domain.getProcessorDomain();
		const std::string log_file_name = (popt.OutputFile + "_ireal_" + num2Padstr(ireal, 4) + "_iter_" + num2Padstr(iter, 4) + "_proc_" + num2Padstr(my_rank, 4) + ".traj");

		std::ofstream log_file;
		log_file.open(log_file_name.c_str());

		int trace_iter = 0;
		std::vector<Streamline> Snew;
		while (true) {
			Snew.clear();
			for (unsigned int i = 0; i < S.size(); ++i) {
				bool tf;
				Domain.bisInProcessorPolygon(S[i].getLastParticle().getP(),tf);
				if (!tf)
					continue;
				//std::cout << "Proc " << world.rank() << " :Particle id: " << S[i].getSid() << std::endl;
				ExitReason er = traceInner(S[i]);
                //std::cout << "Proc " << world.rank() << " :Particle id: " << S[i].getSid() << " " << castExitReasons2String(er) << std::endl;

				if (er == ExitReason::FIRST_POINT_GHOST || er == ExitReason::FAR_AWAY){
                    WRITE::PrintExitReason(log_file, S[i], er);
                    continue;
				}

				for (unsigned int j = 0; j < S[i].size()-1; ++j) {
					WRITE::PrintParticle2Log(log_file, S[i], j);
				}

				if (er == ExitReason::CHANGE_PROCESSOR) {
					Snew.push_back(Streamline(
						S[i].getEid(), S[i].getSid(), S[i].getLastParticle(), S[i].getBBlow(), S[i].getBBupp(), S[i].StuckIter(), S[i].getAge()));
				}
				else {
					WRITE::PrintParticle2Log(log_file, S[i], S[i].size() - 1);
					WRITE::PrintExitReason(log_file, S[i], er);
				}
			}

			world.barrier();
			//std::cout << "Proc " << my_rank << " will send " << Snew.size() << " particles" << std::endl;
			int N_part2track = Snew.size();
			//std::cout << "Proc " << my_rank << " N_part2send = " << N_part2send << std::endl;
			MPI::sumScalar<int>(N_part2track, n_proc, world, MPI_INT);
			//std::cout << "Proc (after) " << my_rank << " N_part2send = " << N_part2send << std::endl;
			world.barrier();
			
			if (my_rank == 0)
				std::cout << "||		Number of active particles: " << N_part2track << " --------" << std::endl << std::flush;

			if (N_part2track == 0)
				break;

			MPI::Send_receive_streamlines(Snew, S, thisDomainPolygon, world);
			//std::cout << my_rank << " S new size: " << S.size() << std::endl;
			
			//if (trace_iter == 10) {
			//	std::cout << "exit after " << trace_iter << " Iteration. DONT FORGET TO REMOVE THIS Condition" << std::endl;
			//	break;
			//}
			trace_iter++;
			//std::cout << trace_iter << std::endl;
			if (trace_iter > popt.MaxProcessorExchanges) {
				if (my_rank == 0) {
					std::cout << "@%$^@&#$^@ exit after MaxProcessorExchanges was reached: trace_iter = " << trace_iter << std::endl;
				}
				break;
			}
		}
		log_file.close();
		//std::cout << "Is here" << std::endl;
	}

	ExitReason ParticleTrace::traceInner(Streamline& S) {
		// Reset the step size
		adaptStepSize = popt.StepOpt.StepSize;
		actualStep = popt.StepOpt.StepSize;
		actualStepTime = popt.StepOpt.StepSizeTime;
		VF.reset();
		XYZ.reset();
		vec3 v;
		vec3 p = S.getLastParticle().getP();
		double tm = S.getLastParticle().getTime();
		int proc = world.rank();
		// Check if the particle is in the domain and exit immediately if not
		ExitReason er = CheckNewPoint(p);
		if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP)
			return ExitReason::INIT_OUT;

		// Attempt to calculate the velocity
        {
            std::vector<int> ids;
            std::vector<double> weights;
            std::map<int, double> proc_map;
            std::vector<vec3> tmp_data;
            bool tf;
            pp = p;
            XYZ.calcWeights(p,ids,weights,proc_map,ll,uu,tf);
            if (!tf) {
                return ExitReason::FAR_AWAY;
            }
            VF.calcVelocity(v,ids,weights,tm);
            //VF.getVec3Data(tmp_data);
        }

		//proc = VF.calcProcID(proc_map);
		//if (proc != world.rank())
		//	return ExitReason::FIRST_POINT_GHOST;

		// If the particle id is 0 then we have to assign velocity
		// If its not zero that means that the velocity has been calculated by the processor
		// that was tracking this streamline before end up in the current processor
		if (S.getLastParticle().getPid() == 0)
			S.InitialVelocity(v);


		//S.getLastParticle().displayAsVEX(true);

		int count_iterations = 0;
		//DEBUG::displayParticleasVex(S.getLastParticle(), true);
		while (er == ExitReason::NO_EXIT) {
			//if (S.size()>1){// Check if there is rapid change in the d
			//	vec3 v_prev = S.getParticleBeforeLast().getV();
			//	vec3 v_curr = S.getLastParticle().getV();
			//	double angle = v_prev.angle(v_curr);
			//}
			if (popt.UpdateStepSize == 1 && popt.method != SolutionMethods::RK45){
			    vv = v;
                updateStep(S.getAge());
			}

			bool foundPoint = findNextPoint(S.getLastParticle(), p, v, tm, er);
			if (foundPoint) {
				//er = CheckNewPointAndCalcVelocity(p, v/*, proc*/, tm);
				//ICHNOS::DBG::displayPVasVex(p, v);
				S.AddParticle(Particle(p, v, S.getLastParticle(), tm));
				//S.getLastParticle().displayAsVEX(true);
				// In the unlike event that the velocity of the point is zero
				// we can avoid unnessecary iterations by simply exit
				if (v.isZero() & (er == ExitReason::NO_EXIT))
					return ExitReason::STUCK;

				bool isNearAttractor = false;
				Domain.bisInNearAttractor(p, isNearAttractor);
				if (isNearAttractor)
					return ExitReason::ATTRACT;
			}

			count_iterations++;
			if (count_iterations > popt.MaxIterationsPerStreamline)
				return ExitReason::MAX_INNER_ITER;
			if (S.StuckIter() > popt.StuckIterations)
				return ExitReason::STUCK;
			if (popt.AgeLimit > 0) {
				if (S.getAge() > popt.AgeLimit)
					return ExitReason::MAX_AGE;
			}
		}
		S.Close(er);
		return er;
	}

	ExitReason ParticleTrace::CheckNewPointAndCalcVelocity(vec3& p, vec3& v/*, int& proc*/, double tm) {
		// This function checks only if the point is still inside the domain
		// It does not check if its in the velocity field owned by this processor
		ExitReason exitreason = CheckNewPoint(p);
		std::map<int, double> proc_map;
		std::vector<int> ids;
		std::vector<double> weights;
        std::vector<vec3> tmp_vec;
		bool tf;
		// If the point is outside the domain we can still calculate velocity if
		// the velocity field allows it. However the particle tracking will be terminated
		if (exitreason != ExitReason::NO_EXIT) {
			if (VF.InterpolateOutsideDomain) {
			    pp = p;
			    XYZ.calcWeights(p,ids,weights,proc_map, ll, uu, tf);
			    //XYZ.sendVec3Data(tmp_vec);
                if (!tf) {
                    v = vec3();
                }
                else{
                    VF.calcVelocity(v,ids,weights, tm);
                    //VF.getVec3Data(tmp_vec);
                }
			}
			else {
				v = vec3();
			}
			return exitreason;
		}
		else {
		    pp = p;
            XYZ.calcWeights(p, ids, weights, proc_map, ll, uu, tf);
            //XYZ.sendVec3Data(tmp_vec);
            if (!tf) {
                v = vec3(-99999, -99999, -99999);
            }
            else{
                VF.calcVelocity(v,ids,weights, tm);
                //VF.getVec3Data(tmp_vec);
            }

			if (v.isInvalid())
				return ExitReason::FAR_AWAY;

			Domain.bisInProcessorPolygon(p, tf);
			if (tf) { // If the point is still in the processor domain then no reason to change even if the other processors may have greater influence
				//proc = world.rank();
				return ExitReason::NO_EXIT;
			}
			// If the point is outside the actual domain then calculate the ownership threshold
			tf = VF.bIsInGhostArea(proc_map);
			if (tf) {// If the point is influenced more by points of another processor then change
				return ExitReason::CHANGE_PROCESSOR;
			}
			else {// If the point is in the overlapping domain but this processor has greater influence
				// make sure that the particle is not outside the expanded domain.
				// If it is change immediately
				Domain.bisInExpandedPolygon(p, tf);
				if (!tf) {
					return ExitReason::CHANGE_PROCESSOR;
				}
				else{
					//proc = world.rank();
					return ExitReason::NO_EXIT;
				}
			}
		}
	}

	ExitReason ParticleTrace::CheckNewPoint(vec3& p) {
		bool bInDomain = true;
		Domain.bisInPolygon(p, bInDomain);
		if (!bInDomain)
			return ExitReason::EXIT_SIDE;

		double top, bottom;
		Domain.getTopBottomElevation(p, top, bottom, bInDomain);
		if (!bInDomain)
            return ExitReason::FAR_AWAY;

		if (p.z < bottom)
			return ExitReason::EXIT_BOTTOM;
		if (p.z > top)
			return ExitReason::EXIT_TOP;

		return ExitReason::NO_EXIT;
	}

	bool ParticleTrace::findNextPoint(const Particle& P, vec3& pnew, vec3& vnew, double& tm, ExitReason& er) {
		switch (popt.method)
		{
		case SolutionMethods::Euler:
			return EulerStep(P, pnew, vnew, tm, er);
			break;
		case SolutionMethods::RK2:
			return RK2Step(P, pnew, vnew, tm, er);
			break;
		case SolutionMethods::RK4:
			return RK4Step(P, pnew, vnew, tm, er);
			break;
		case SolutionMethods::RK45:
			return RK45Step(P, pnew, vnew, tm, er);
			break;
		default:
			return EulerStep(P, pnew, vnew, tm, er);
			break;
		}
	}

	bool ParticleTrace::EulerStep(const Particle& P, vec3& pnew, vec3& vnew, double& tm, ExitReason& er) {
	    double init_tm = tm;
		//DEBUG::displayParticleasVex(P, true);
		pnew = P.getP() + P.getV().normalize() * actualStep * popt.Direction;
		//DEBUG::displayVectorasVex(pnew);
        tm += actualStep*popt.Direction / P.getV().len();

        er = CheckNewPointAndCalcVelocity(pnew, vnew/*, proc*/, tm);
        if (actualStep > popt.StepOpt.minExitStepSize){
            if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP ){
                actualStep = actualStep * 0.45;
                tm = init_tm;
                return EulerStep(P, pnew, vnew, tm, er);
            }
        }
		return true;
	}

	bool ParticleTrace::RK2Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, ExitReason& er) {
        double init_tm = tm;
		//DEBUG::displayParticleasVex(P, true);
		int proc = world.rank();
		// Find the new position using the velocity of the current point
		vec3 p1 = P.getV().normalize() * actualStep * popt.Direction + P.getP();
        double tm1 = tm + actualStep*popt.Direction / P.getV().len();
		
		// Find the velocity on that point
		vec3 v1;
		er = CheckNewPointAndCalcVelocity(p1, v1/*, proc*/, tm1);
		//DEBUG::displayPVasVex(p1, v1);
		if (v1.isZero()) return false;
		// Take a step with the average velocity
		vec3 v_m = v1 * 0.5 + P.getV() * 0.5;
		pnew = P.getP() + v_m.normalize() * actualStep* popt.Direction;
		//DEBUG::displayVectorasVex(pnew);
        tm += actualStep*popt.Direction / v_m.len();
        er = CheckNewPointAndCalcVelocity(pnew, vnew/*, proc*/, tm);
        if (actualStep > popt.StepOpt.minExitStepSize){
            if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP ){
                actualStep = actualStep * 0.45;
                tm = init_tm;
                return RK2Step(P, pnew, vnew, tm, er);
            }
        }
		return true;
	}

	bool ParticleTrace::RK4Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, ExitReason& er) {
        double init_tm = tm;
		//DEBUG::displayParitcleasVex(P, true);
		// See RK2Step
		int proc = world.rank();
		// Step 1 ----------------------------
		// take a half step from the current point using the velocity of that point
		vec3 p2 = P.getV().normalize() * (actualStep/2) * popt.Direction + P.getP();
        double tm2 = tm + (actualStep/2)*popt.Direction / P.getV().len();
		// Find the velocity on point p2 that point
		vec3 v2;
		er = CheckNewPointAndCalcVelocity(p2, v2/*, proc*/, tm2);
		//DEBUG::displayPVasVex(p2, v2);
		if (v2.isZero()) return false;

		// Step 2 ----------------------------
		// take a half step from the current point using the v2 velocity
		vec3 p3 = v2.normalize() * (actualStep/2) * popt.Direction + P.getP();
        double tm3 = tm + (actualStep/2)*popt.Direction / v2.len();
		// Find the velocity on point p3 that point
		vec3 v3;
		er = CheckNewPointAndCalcVelocity(p3, v3/*, proc*/, tm3);
		//DEBUG::displayPVasVex(p3, v3);
		if (v3.isZero()) return false;

		// Step 3 ----------------------------
		// take a full step from the current point using the v3 velocity
		vec3 p4 = v3.normalize() * actualStep * popt.Direction + P.getP();
        double tm4 = tm + actualStep*popt.Direction / v3.len();
		// Find the velocity on point p4 that point
		vec3 v4;
		er = CheckNewPointAndCalcVelocity(p4, v4/*, proc*/, tm4);
		//DEBUG::displayPVasVex(p4, v4);
		if (v4.isZero()) return false;

		vec3 v_m = (P.getV() + v2 * 2.0 + v3 * 2.0 + v4)*(1.0/6.0);
		pnew = v_m.normalize() * actualStep * popt.Direction + P.getP();
		//DEBUG::displayVectorasVex(pnew);
        tm += actualStep*popt.Direction / v_m.len();
        er = CheckNewPointAndCalcVelocity(pnew, vnew/*, proc*/, tm);
        if (actualStep > popt.StepOpt.minExitStepSize){
            if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP ){
                actualStep = actualStep * 0.45;
                tm = init_tm;
                return RK4Step(P, pnew, vnew, tm, er);
            }
        }
		return true;
	}

	bool ParticleTrace::RK45Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, ExitReason& er) {
        double init_tm = tm;
		//DEBUG::displayParticleasVex(P, true);
		// See RK2Step
		int proc = world.rank();
		vec3 v_m, p2, p3, p4, p5, p6, v2, v3, v4, v5, v6, yn, zn;
		vec3 v1 = P.getV();

		double ssdir = adaptStepSize * popt.Direction;

		//Runge-Kutta-Fehlberg Method (RKF45)
		// Step 1 ----------------------------
		int istep = 0;
		p2 = v1.normalize() * CF[istep][0] * ssdir + P.getP();
		double tm2 = tm + CF[istep][0] * ssdir * v1.len();
		er = CheckNewPointAndCalcVelocity(p2, v2/*, proc*/, tm2);
		//DEBUG::displayPVasVex(p2, v2);
		if (v2.isZero()) { return false; }
		if (adaptStepSize > popt.StepOpt.minExitStepSize)
			if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP) {
				adaptStepSize = adaptStepSize * 0.24;
                tm = init_tm;
				return RK45Step(P, pnew, vnew, tm, er);
			}

		// Step 2 ----------------------------
		istep++;
		v_m = v1 * CF[istep][1] + v2 * CF[istep][2];
		p3 = v_m.normalize() * CF[istep][0] * ssdir + P.getP();
		double tm3 = tm + CF[istep][0] * ssdir * v_m.len();
		er = CheckNewPointAndCalcVelocity(p3, v3/*, proc*/, tm3);
		//DEBUG::displayPVasVex(p3, v3);
		if (v3.isZero()) return false;
		if (adaptStepSize > popt.StepOpt.minExitStepSize)
			if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP) {
				adaptStepSize = adaptStepSize * 0.37;
                tm = init_tm;
                return RK45Step(P, pnew, vnew, tm, er);
			}

		// Step 3 ----------------------------
		istep++;
		v_m = v1 * CF[istep][1] + v2 * CF[istep][2] + v3 * CF[istep][3];
		p4 = v_m.normalize() * CF[istep][0] * ssdir + P.getP();
		double tm4 = tm + CF[istep][0] * ssdir * v_m.len();
		er = CheckNewPointAndCalcVelocity(p4, v4/*, proc*/, tm4);
		//DEBUG::displayPVasVex(p4, v4);
		if (v4.isZero()) return false;

		// Step 4 ----------------------------
		istep++;
		v_m = v1 * CF[istep][1] + v2 * CF[istep][2] + v3 * CF[istep][3] + v4 * CF[istep][4];
		p5 = v_m.normalize() * CF[istep][0] * ssdir + P.getP();
        double tm5 = tm + CF[istep][0] * ssdir * v_m.len();
		er = CheckNewPointAndCalcVelocity(p5, v5/*, proc*/, P.getTime());
		//DEBUG::displayPVasVex(p5, v5);
		if (v5.isZero()) return false;

		// Step 5 ----------------------------
		istep++;
		v_m = v1 * CF[istep][1] + v2 * CF[istep][2] + v3 * CF[istep][3] + v4 * CF[istep][4] + v5 * CF[istep][5];
		p6 = v_m.normalize() * CF[istep][0] * ssdir + P.getP();
        double tm6 = tm + CF[istep][0] * ssdir * v_m.len();
		er = CheckNewPointAndCalcVelocity(p6, v6/*, proc*/, tm6);
		//DEBUG::displayPVasVex(p6, v6);
		if (v6.isZero()) return false;
		if ( adaptStepSize > popt.StepOpt.minExitStepSize )
			if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP ) {
				adaptStepSize = adaptStepSize * 0.45;
                tm = init_tm;
                return RK45Step(P, pnew, vnew, tm, er);
			}

		//// Calculate the point using two different paths
		istep++;
		v_m = v1 * CF[istep][1] + v3 * CF[istep][2] + v4 * CF[istep][3] + v5 * CF[istep][4];
		yn = v_m.normalize() * CF[istep][0] * ssdir + P.getP();
		//DEBUG::displayVectorasVex(yn);
		istep++;
		v_m = v1 * CF[istep][1] + v3 * CF[istep][2] + v4 * CF[istep][3] + v5 * CF[istep][4] + v6 * CF[istep][5];
		zn = v_m.normalize() * CF[istep][0] * ssdir + P.getP();
		//DEBUG::displayVectorasVex(zn);
		double R = (yn - zn).len();
		double q = 0.84 * pow(popt.AdaptOpt.ToleranceStepSize / R, 0.25);

		if (q >= 1) {
			q = std::min(q, popt.AdaptOpt.increaseRateChange);
			// We can increase the step size
			adaptStepSize = adaptStepSize * q;
			if (adaptStepSize > popt.AdaptOpt.MaxStepSize)
				adaptStepSize = popt.AdaptOpt.MaxStepSize;
			//if (adaptStepSize < popt.MinStepSize)
			//	adaptStepSize = popt.MinStepSize;
			pnew = yn;
            tm += CF[istep][0] * ssdir / v_m.len();
		}
		else {
			if (adaptStepSize <= popt.AdaptOpt.MinStepSize) {
				pnew = yn;
                tm += CF[istep][0] * ssdir / v_m.len();
				if (R > popt.AdaptOpt.ToleranceStepSize){
					std::cout << "The algorithm will continue although it cannot reach the desired accuracy of "
						<< popt.AdaptOpt.ToleranceStepSize
						<< "with the limit of minimum step size " << popt.AdaptOpt.MinStepSize << std::endl;
					std::cout << "Either relax the tolerance or reduce the Minimum Step Size" << std::endl;
				}
			}
			else {
				q = std::min(q, popt.AdaptOpt.limitUpperDecreaseStep);
				adaptStepSize = adaptStepSize * q;
				if (R > popt.AdaptOpt.ToleranceStepSize){
                    tm = init_tm;
                    return RK45Step(P, pnew, vnew, tm, er);
				}
			}
		}
        er = CheckNewPointAndCalcVelocity(pnew, vnew/*, proc*/, tm);
        if (adaptStepSize > popt.StepOpt.minExitStepSize){
            if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP ){
                adaptStepSize = adaptStepSize * 0.45;
                tm = init_tm;
                return RK45Step(P, pnew, vnew, tm, er);
            }
        }
		return true;
	}

	void ParticleTrace::updateStep(double currentAge) {
        double stepLen, stepTime;
        double dst = ICHNOS::diameter_along_velocity(pp, vv, ll, uu);

        if (popt.StepOpt.nSteps > 0){
            // This is the length that respects the nSteps
            stepLen = dst/popt.StepOpt.nSteps;
        }
        else{
            stepLen = 10000000;
        }

        if (popt.StepOpt.StepSize > 0){
            // If the length of the user step is smaller than the step defined by the nSteps.
            // Use the used defined that is smaller
            if (popt.StepOpt.StepSize < stepLen){
                stepLen = popt.StepOpt.StepSize;
            }
        }

        double vlen = vv.len();
        if (popt.StepOpt.StepSizeTime > 0){
            stepTime = vlen*popt.StepOpt.StepSizeTime;
        }
        else{
            stepTime = 10000000;
        }

        if (popt.AgeLimit > 0){
            double rem_age = popt.AgeLimit - currentAge;
            rem_age = vlen*(rem_age + 0.1);
            stepTime = std::min<double>(stepTime, rem_age);
        }
        //TODO Add the part of the nStepsTime limit
        actualStep = std::min<double>(stepTime, stepLen);
	}
}
