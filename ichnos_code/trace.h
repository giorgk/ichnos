#pragma once
#include <boost/mpi.hpp>
#include <boost/multi_array.hpp>

#include <thread>

#if _USEHF > 0
#include <highfive/H5File.hpp>
#endif

#include "ichnos_structures.h"
#include "velocity_base.h"
#include "ichnos_domain.h"
#include "ichnos_XYZ_base.h"


namespace ICHNOS {

	class ParticleTrace {
	public:
		ParticleTrace(
			boost::mpi::communicator& world_in,
            //XYZ_base &XYZ_in,
			velocityField &VF_in,
			DomainBase &dom_in,
			ParticleOptions& popt_in);

		void Trace();
		void Trace(int ireal);
		void multi_threading_trace();
		void run_thread(int id);

		bool readInputFiles(std::vector<Streamline>& S);

	private:
		const ParticleOptions popt;
		velocityField& VF;
		DomainBase& Domain;
        //XYZ_base& XYZ;
		boost::mpi::communicator world;
		bool findNextPoint(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu,  ExitReason& er);
		ExitReason traceInner(Streamline& S);
		void traceOuter(std::vector<Streamline>& S, int iter, int ireal);
		ExitReason CheckNewPoint(vec3& p);
		ExitReason CheckNewPointAndCalcVelocity(vec3& p, vec3& v,helpVars& pvlu/*, int &proc*/, double tm = 0);

		bool EulerStep(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu, ExitReason& er);
		bool RK2Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu, ExitReason& er);
		bool RK4Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu, ExitReason& er);
		bool RK45Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu, ExitReason& er);
        bool RalstonStep(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu, ExitReason& er);
        bool PECEStep(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu, ExitReason& er);

		void updateStep(double currentAge, helpVars& pvlu);

		const std::vector<std::vector<double> > CF = getRK45coef();

		//double adaptStepSize;
		//double actualStep;
		//double actualStepTime;
		//vec3 pp, vv, ll, uu;

        std::vector<Streamline> Streamlines4thread;

	};

	ParticleTrace::ParticleTrace(
		boost::mpi::communicator& world_in,
        //XYZ_base &XYZ_in,
		velocityField& VF_in,
		DomainBase& dom_in,
		ParticleOptions &popt_in)
		:
		world(world_in),
        //XYZ(XYZ_in),
		VF(VF_in),
		Domain(dom_in),
		popt(popt_in)
	{
		//adaptStepSize = popt.StepOpt.StepSize;
	}

	bool ParticleTrace::readInputFiles(std::vector <Streamline>& S) {
		bool tfw = false;
		bool tfp = false;
		if (!popt.WellFile.empty()) {
			tfw = READ::readWellFile(popt.WellFile, S, VF.isVelTransient());
		}

		if (!popt.ParticleFile.empty()) {
			tfp = READ::readParticleFile(popt.ParticleFile, S, VF.isVelTransient());
		}

		if (tfw || tfp)
			return true;
		else
			return false;
	}

	void ParticleTrace::Trace() {
	    if (world.size() == 1 || popt.RunAsThread){
	        multi_threading_trace();
	    }
	    else{
            for (int ireal = 0; ireal < popt.Nrealizations; ++ireal){
                Trace(ireal);
            }
	    }
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
		//const std::string log_file_name = (popt.OutputFile + "_ireal_" + num2Padstr(ireal, 4) + "_iter_" + num2Padstr(iter, 4) + "_proc_" + num2Padstr(my_rank, 4) + ".traj");

		//std::ofstream log_file;
		//log_file.open(log_file_name.c_str());
        bool append = false;
		int trace_iter = 0;
		std::vector<Streamline> Snew;
		while (true) {
			Snew.clear();
			for (unsigned int i = 0; i < S.size(); ++i) {
				bool tf;
				Domain.bisInProcessorPolygon(S[i].getLastParticle().getP(),tf);
				if (!tf){
                    S[i].Close(ExitReason::NOT_IN_SUBDOMAIN);
                    continue;
                }

				//std::cout << "Proc " << world.rank() << " :Particle id: " << S[i].getSid() << std::endl;
				ExitReason er = traceInner(S[i]);
				S[i].Close(er);
                //std::cout << "Proc " << world.rank() << " :Particle id: " << S[i].getSid() << " " << castExitReasons2String(er) << std::endl;

				if (er == ExitReason::CHANGE_PROCESSOR) {
					Snew.push_back(Streamline(
						S[i].getEid(), S[i].getSid(), S[i].getLastParticle(), S[i].getBBlow(), S[i].getBBupp(), S[i].StuckIter(), S[i].getAge()));
				}

			}
            const std::string out_file_name = (popt.OutputFile + "_ireal_" + num2Padstr(ireal, 4) + "_iter_" + num2Padstr(iter, 4) + "_proc_" + num2Padstr(my_rank, 4));
            WRITE::writeStreamlines(S, out_file_name, popt.printH5, popt.printASCII, append);
            append = true;

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
		//log_file.close();
		//std::cout << "Is here" << std::endl;
	}

	ExitReason ParticleTrace::traceInner(Streamline& S) {
	    //if (S.getSid() == 2){
	    //    std::cout << "Hi" << std::endl;
	    //}
		// Reset the step size
		S.PVLU.adaptStepSize = popt.StepOpt.StepSize;
        S.PVLU.actualStep = popt.StepOpt.StepSize;
        S.PVLU.actualStepTime = popt.StepOpt.StepSizeTime;
		VF.reset(S);
		VF.XYZ.reset(S);
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
            std::map<int, double> proc_map;
            bool tf;
            S.PVLU.pp = p;

            VF.calcVelocity(p, v, proc_map, S.PVLU, tf, tm);
            if (!tf){
                return ExitReason::FIRST_POINT_GHOST;
            }
            //std::cout << world.rank() << ": " << v.x << "," << v.y << "," << v.z << std::endl;
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
			    S.PVLU.vv = v;
                updateStep(S.getAge(), S.PVLU);
			}

			bool foundPoint = findNextPoint(S.getLastParticle(), p, v, tm, S.PVLU, er);
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
		return er;
	}

	ExitReason ParticleTrace::CheckNewPointAndCalcVelocity(vec3& p, vec3& v, helpVars& pvlu/*, int& proc*/, double tm) {
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
                pvlu.pp = p;
			    //VF.XYZ.calcWeights(p,ids,weights,proc_map, pvlu, tf);
			    //XYZ.sendVec3Data(tmp_vec);
                VF.calcVelocity(p, v, proc_map, pvlu, tf, tm);
                if (!tf) {
                    v = vec3();
                }
                //else{

                    //VF.getVec3Data(tmp_vec);
                //}
			}
			else {
				v = vec3();
			}
			return exitreason;
		}
		else {
            pvlu.pp = p;
            //VF.XYZ.calcWeights(p, ids, weights, proc_map, pvlu, tf);
            //XYZ.sendVec3Data(tmp_vec);
            VF.calcVelocity(p, v, proc_map, pvlu, tf, tm);
            if (!tf) {
                v = vec3(-99999, -99999, -99999);
            }
            //else{

                //VF.getVec3Data(tmp_vec);
            //}

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
				//Domain.bisInExpandedPolygon(p, tf);
				//if (!tf) {
				//	return ExitReason::CHANGE_PROCESSOR;
				//}
				//else{
					//proc = world.rank();
					return ExitReason::NO_EXIT;
				//}
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
        //std::cout << p.z << ":" << top << "-" << bottom << std::endl;
		if (!bInDomain)
            return ExitReason::FAR_AWAY;

		if (p.z < bottom)
			return ExitReason::EXIT_BOTTOM;
		if (p.z > top)
			return ExitReason::EXIT_TOP;

		return ExitReason::NO_EXIT;
	}

	bool ParticleTrace::findNextPoint(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu, ExitReason& er) {
		switch (popt.method)
		{
		case SolutionMethods::Euler:
			return EulerStep(P, pnew, vnew, tm, pvlu, er);
			break;
		case SolutionMethods::RK2:
			return RK2Step(P, pnew, vnew, tm, pvlu, er);
			break;
		case SolutionMethods::RK4:
			return RK4Step(P, pnew, vnew, tm, pvlu, er);
			break;
		case SolutionMethods::RK45:
			return RK45Step(P, pnew, vnew, tm, pvlu, er);
			break;
        case SolutionMethods::RAL:
            return RalstonStep(P, pnew, vnew, tm, pvlu, er);
            break;
        case SolutionMethods::PECE:
            return PECEStep(P, pnew, vnew, tm, pvlu, er);
            break;
		default:
			return EulerStep(P, pnew, vnew, tm, pvlu, er);
			break;
		}
	}

	bool ParticleTrace::EulerStep(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu, ExitReason& er) {
	    double init_tm = tm;
		//DEBUG::displayParticleasVex(P, true);
		pnew = P.getP() + P.getV().normalize() * pvlu.actualStep * popt.Direction;
		//DEBUG::displayVectorasVex(pnew);
        tm += pvlu.actualStep*popt.Direction / P.getV().len();

        er = CheckNewPointAndCalcVelocity(pnew, vnew, pvlu/*, proc*/, tm);
        if (pvlu.actualStep > popt.StepOpt.minExitStepSize && popt.UpdateStepSize == 1){
            if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP ){
                pvlu.actualStep = pvlu.actualStep * 0.45;
                tm = init_tm;
                return EulerStep(P, pnew, vnew, tm, pvlu, er);
            }
        }
		return true;
	}

    bool ParticleTrace::PECEStep(const Particle &P, vec3 &pnew, vec3 &vnew, double &tm, helpVars &pvlu,
                                 ExitReason &er) {
        bool tf;
        double init_tm = tm;
        double tmp_tm = tm;
        // Take an Euler step
        vec3 vn_1, vn, pn_1, pn;
        tf = EulerStep(P,pn, vn, tmp_tm, pvlu, er);
        if (!tf) return false;
        for (int i = 0; i < popt.PECEOpt.Order; ++i){

            vn_1 = P.getV()*0.5 + vn*0.5;

            pn_1 =  P.getP() + vn_1.normalize() * pvlu.actualStep * popt.Direction;
            if (pn_1.distance(pn.x, pn.y, pn.z) < popt.PECEOpt.Tolerance){
                break;
            }
            else{
                pn = pn_1;
                vn = vn_1;
            }
        }
    }

    bool ParticleTrace::RalstonStep(const Particle &P, vec3 &pnew, vec3 &vnew, double &tm, helpVars &pvlu, ExitReason &er) {
        double init_tm = tm;

        // Take a 2/3 step from the current point using the velocity of the current point
        vec3 p2 = P.getV().normalize() * (2.0*pvlu.actualStep/3.0) * popt.Direction + P.getP();
        double tm2 = tm + (2.0*pvlu.actualStep/3.0)*popt.Direction / P.getV().len();
        // Find the velocity on point p2 that point
        vec3 v2;
        er = CheckNewPointAndCalcVelocity(p2, v2, pvlu/*, proc*/, tm2);
        //DEBUG::displayPVasVex(p2, v2);
        if (v2.isZero()) return false;

        // take another 2/3 step from the current point using the v2 velocity
        vec3 p3 = v2.normalize() * (2.0*pvlu.actualStep/3.0) * popt.Direction + P.getP();
        double tm3 = tm + (2.0*pvlu.actualStep/3.0)*popt.Direction / v2.len();
        // Find the velocity on point p3 that point
        vec3 v3;
        er = CheckNewPointAndCalcVelocity(p3, v3, pvlu/*, proc*/, tm3);
        //DEBUG::displayPVasVex(p3, v3);
        if (v3.isZero()) return false;

        vec3 v_m = P.getV()*0.25 + v3*0.75 ;

        pnew = v_m.normalize() * pvlu.actualStep * popt.Direction + P.getP();
        //DEBUG::displayVectorasVex(pnew);
        tm += pvlu.actualStep*popt.Direction / v_m.len();
        er = CheckNewPointAndCalcVelocity(pnew, vnew, pvlu/*, proc*/, tm);
        if (pvlu.actualStep > popt.StepOpt.minExitStepSize){
            if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP ){
                pvlu.actualStep = pvlu.actualStep * 0.45;
                tm = init_tm;
                return RalstonStep(P, pnew, vnew, tm, pvlu, er);
            }
        }
        return true;
    }

	bool ParticleTrace::RK2Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu, ExitReason& er) {
        double init_tm = tm;
		//DEBUG::displayParticleasVex(P, true);
		int proc = world.rank();
		// Find the new position using the velocity of the current point
		vec3 p1 = P.getV().normalize() * pvlu.actualStep * popt.Direction + P.getP();
        double tm1 = tm + pvlu.actualStep*popt.Direction / P.getV().len();
		
		// Find the velocity on that point
		vec3 v1;
		er = CheckNewPointAndCalcVelocity(p1, v1, pvlu/*, proc*/, tm1);
		//DEBUG::displayPVasVex(p1, v1);
		if (v1.isZero()) return false;
		// Take a step with the average velocity
		vec3 v_m = v1 * 0.5 + P.getV() * 0.5;
		pnew = P.getP() + v_m.normalize() * pvlu.actualStep* popt.Direction;
		//DEBUG::displayVectorasVex(pnew);
        tm += pvlu.actualStep*popt.Direction / v_m.len();
        er = CheckNewPointAndCalcVelocity(pnew, vnew, pvlu/*, proc*/, tm);
        if (pvlu.actualStep > popt.StepOpt.minExitStepSize){
            if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP ){
                pvlu.actualStep = pvlu.actualStep * 0.45;
                tm = init_tm;
                return RK2Step(P, pnew, vnew, tm, pvlu, er);
            }
        }
		return true;
	}

	bool ParticleTrace::RK4Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu, ExitReason& er) {
        double init_tm = tm;
		//DEBUG::displayParitcleasVex(P, true);
		// See RK2Step
		int proc = world.rank();
		// Step 1 ----------------------------
		// take a half step from the current point using the velocity of the current point
		vec3 p2 = P.getV().normalize() * (pvlu.actualStep/2.0) * popt.Direction + P.getP();
        double tm2 = tm + (pvlu.actualStep/2.0)*popt.Direction / P.getV().len();
		// Find the velocity on point p2 that point
		vec3 v2;
		er = CheckNewPointAndCalcVelocity(p2, v2, pvlu/*, proc*/, tm2);
		//DEBUG::displayPVasVex(p2, v2);
		if (v2.isZero()) return false;

		// Step 2 ----------------------------
		// take a half step from the current point using the v2 velocity
		vec3 p3 = v2.normalize() * (pvlu.actualStep/2.0) * popt.Direction + P.getP();
        double tm3 = tm + (pvlu.actualStep/2.0)*popt.Direction / v2.len();
		// Find the velocity on point p3 that point
		vec3 v3;
		er = CheckNewPointAndCalcVelocity(p3, v3, pvlu/*, proc*/, tm3);
		//DEBUG::displayPVasVex(p3, v3);
		if (v3.isZero()) return false;

		// Step 3 ----------------------------
		// take a full step from the current point using the v3 velocity
		vec3 p4 = v3.normalize() * pvlu.actualStep * popt.Direction + P.getP();
        double tm4 = tm + pvlu.actualStep*popt.Direction / v3.len();
		// Find the velocity on point p4 that point
		vec3 v4;
		er = CheckNewPointAndCalcVelocity(p4, v4, pvlu/*, proc*/, tm4);
		//DEBUG::displayPVasVex(p4, v4);
		if (v4.isZero()) return false;

		vec3 v_m = (P.getV() + v2 * 2.0 + v3 * 2.0 + v4)*(1.0/6.0);
		pnew = v_m.normalize() * pvlu.actualStep * popt.Direction + P.getP();
		//DEBUG::displayVectorasVex(pnew);
        tm += pvlu.actualStep*popt.Direction / v_m.len();
        er = CheckNewPointAndCalcVelocity(pnew, vnew, pvlu/*, proc*/, tm);
        if (pvlu.actualStep > popt.StepOpt.minExitStepSize){
            if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP ){
                pvlu.actualStep = pvlu.actualStep * 0.45;
                tm = init_tm;
                return RK4Step(P, pnew, vnew, tm, pvlu, er);
            }
        }
		return true;
	}

	bool ParticleTrace::RK45Step(const Particle& P, vec3& pnew, vec3& vnew, double& tm, helpVars& pvlu, ExitReason& er) {
        double init_tm = tm;
		//DEBUG::displayParticleasVex(P, true);
		// See RK2Step
		int proc = world.rank();
		vec3 v_m, p2, p3, p4, p5, p6, v2, v3, v4, v5, v6, yn, zn;
		vec3 v1 = P.getV();

		double ssdir = pvlu.adaptStepSize * popt.Direction;

		//Runge-Kutta-Fehlberg Method (RKF45)
		// Step 1 ----------------------------
		int istep = 0;
		p2 = v1.normalize() * CF[istep][0] * ssdir + P.getP();
		double tm2 = tm + CF[istep][0] * ssdir * v1.len();
		er = CheckNewPointAndCalcVelocity(p2, v2, pvlu/*, proc*/, tm2);
		//DEBUG::displayPVasVex(p2, v2);
		if (v2.isZero()) { return false; }
		if (pvlu.adaptStepSize > popt.StepOpt.minExitStepSize)
			if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP) {
                pvlu.adaptStepSize = pvlu.adaptStepSize * 0.24;
                tm = init_tm;
				return RK45Step(P, pnew, vnew, tm, pvlu, er);
			}

		// Step 2 ----------------------------
		istep++;
		v_m = v1 * CF[istep][1] + v2 * CF[istep][2];
		p3 = v_m.normalize() * CF[istep][0] * ssdir + P.getP();
		double tm3 = tm + CF[istep][0] * ssdir * v_m.len();
		er = CheckNewPointAndCalcVelocity(p3, v3, pvlu/*, proc*/, tm3);
		//DEBUG::displayPVasVex(p3, v3);
		if (v3.isZero()) return false;
		if (pvlu.adaptStepSize > popt.StepOpt.minExitStepSize)
			if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP) {
                pvlu.adaptStepSize = pvlu.adaptStepSize * 0.37;
                tm = init_tm;
                return RK45Step(P, pnew, vnew, tm, pvlu, er);
			}

		// Step 3 ----------------------------
		istep++;
		v_m = v1 * CF[istep][1] + v2 * CF[istep][2] + v3 * CF[istep][3];
		p4 = v_m.normalize() * CF[istep][0] * ssdir + P.getP();
		double tm4 = tm + CF[istep][0] * ssdir * v_m.len();
		er = CheckNewPointAndCalcVelocity(p4, v4, pvlu/*, proc*/, tm4);
		//DEBUG::displayPVasVex(p4, v4);
		if (v4.isZero()) return false;

		// Step 4 ----------------------------
		istep++;
		v_m = v1 * CF[istep][1] + v2 * CF[istep][2] + v3 * CF[istep][3] + v4 * CF[istep][4];
		p5 = v_m.normalize() * CF[istep][0] * ssdir + P.getP();
        double tm5 = tm + CF[istep][0] * ssdir * v_m.len();
		er = CheckNewPointAndCalcVelocity(p5, v5, pvlu/*, proc*/, P.getTime());
		//DEBUG::displayPVasVex(p5, v5);
		if (v5.isZero()) return false;

		// Step 5 ----------------------------
		istep++;
		v_m = v1 * CF[istep][1] + v2 * CF[istep][2] + v3 * CF[istep][3] + v4 * CF[istep][4] + v5 * CF[istep][5];
		p6 = v_m.normalize() * CF[istep][0] * ssdir + P.getP();
        double tm6 = tm + CF[istep][0] * ssdir * v_m.len();
		er = CheckNewPointAndCalcVelocity(p6, v6, pvlu/*, proc*/, tm6);
		//DEBUG::displayPVasVex(p6, v6);
		if (v6.isZero()) return false;
		if ( pvlu.adaptStepSize > popt.StepOpt.minExitStepSize )
			if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP ) {
                pvlu.adaptStepSize = pvlu.adaptStepSize * 0.45;
                tm = init_tm;
                return RK45Step(P, pnew, vnew, tm, pvlu, er);
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
            pvlu.adaptStepSize = pvlu.adaptStepSize * q;
			if (pvlu.adaptStepSize > popt.AdaptOpt.MaxStepSize)
                pvlu.adaptStepSize = popt.AdaptOpt.MaxStepSize;
			//if (adaptStepSize < popt.MinStepSize)
			//	adaptStepSize = popt.MinStepSize;
			pnew = yn;
            tm += CF[istep][0] * ssdir / v_m.len();
		}
		else {
			if (pvlu.adaptStepSize <= popt.AdaptOpt.MinStepSize) {
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
                pvlu.adaptStepSize = pvlu.adaptStepSize * q;
				if (R > popt.AdaptOpt.ToleranceStepSize){
                    tm = init_tm;
                    return RK45Step(P, pnew, vnew, tm, pvlu, er);
				}
			}
		}
        er = CheckNewPointAndCalcVelocity(pnew, vnew, pvlu/*, proc*/, tm);
        if (pvlu.adaptStepSize > popt.StepOpt.minExitStepSize){
            if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP ){
                pvlu.adaptStepSize = pvlu.adaptStepSize * 0.45;
                tm = init_tm;
                return RK45Step(P, pnew, vnew, tm, pvlu, er);
            }
        }
		return true;
	}

	void ParticleTrace::updateStep(double currentAge, helpVars& pvlu) {
        double stepLen, stepTime;
        double dst = ICHNOS::diameter_along_velocity(pvlu.pp, pvlu.vv, pvlu.ll, pvlu.uu);

        if (popt.StepOpt.nSteps >= 1.0){
            // This is the length that respects the nSteps
            stepLen = dst/popt.StepOpt.nSteps;
        }
        else{
            stepLen = 10000000;
        }

        if (popt.StepOpt.StepSize > 0){
            // If the length of the user step is smaller than the step calculated by the nSteps.
            // Use the used defined that is smaller
            if (popt.StepOpt.StepSize < stepLen){
                stepLen = popt.StepOpt.StepSize;
            }
        }

        double vlen = pvlu.vv.len();
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

        if (popt.StepOpt.nStepsTime > 0 && VF.isVelTransient()){
            double stepTimetrans = VF.stepTimeupdate(pvlu);
            if (stepTimetrans < stepTime){
                stepTime = stepTimetrans;
            }
        }
        pvlu.actualStep = std::min<double>(stepTime, stepLen);
	}

	void ParticleTrace::multi_threading_trace() {
        std::vector<Streamline> AllStreamlines;
        bool tf = readInputFiles(AllStreamlines);
        if (!tf){
            return;
        }

        int particle_index_start = 0;
        int particle_index_end = 0;
        int iter = 0;
        double Simulation_time = 0;
        while (true){
            Streamlines4thread.clear();
            particle_index_end = particle_index_start + popt.ParticlesInParallel;
            if (particle_index_end + static_cast<int>(static_cast<double>(popt.ParticlesInParallel)*0.3) > static_cast<int>(AllStreamlines.size())){
                particle_index_end = AllStreamlines.size();
            }

            for (int ireal = 0; ireal < popt.Nrealizations; ++ireal){
                auto start = std::chrono::high_resolution_clock::now();
                if (world.size() > 1 || (world.size() == 1 && popt.Nthreads == 1)){
                    for (int i = particle_index_start + world.rank(); i < particle_index_end; i = i + world.size()){
                        Streamlines4thread.push_back(AllStreamlines[i]);
                    }
                    run_thread(world.rank());

                    auto finish = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = finish - start;
                    Simulation_time = Simulation_time + elapsed.count();

                    std::string filename;
                    if (popt.Nrealizations > 1)
                        filename = popt.OutputFile + "_iter_" + num2Padstr(iter, 4) + "_ireal_" + num2Padstr(ireal, 4) + "_iproc_" + num2Padstr(world.rank(), 4);
                    else{
                        if (world.size() > 1){
                            filename = popt.OutputFile+ "_iter_" + num2Padstr(iter, 4) + "_iproc_" + num2Padstr(world.rank(), 4);
                        }
                        else if(world.size() == 1){
                            filename = popt.OutputFile+ "_iter_" + num2Padstr(iter, 4);
                        }
                    }

                    WRITE::writeStreamlines(Streamlines4thread, filename, popt.printH5, popt.printASCII);
                    world.barrier();
                }
                else{
                    for (int i = particle_index_start; i < particle_index_end; ++i){
                        Streamlines4thread.push_back(AllStreamlines[i]);
                    }
                    std::vector<std::thread> T;
                    for (int i = 0; i < popt.Nthreads; ++i){
                        T.push_back(std::thread(&ParticleTrace::run_thread, this, i));
                    }

                    for (int i = 0; i < popt.Nthreads; ++i){
                        if (T[i].joinable()){
                            T[i].join();
                            //std::cout << "Join" << i << std::endl;
                        }
                    }
                    auto finish = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = finish - start;
                    Simulation_time = Simulation_time + elapsed.count();
                    std::cout << "Printing files" << std::endl;
                    std::string filename;
                    if (popt.Nrealizations > 1)
                        filename = popt.OutputFile + "_iter_" + num2Padstr(iter, 4) + "_ireal_" + num2Padstr(ireal, 4);
                    else
                        filename = popt.OutputFile+ "_iter_" + num2Padstr(iter, 4);
                    WRITE::writeStreamlines(Streamlines4thread, filename, popt.printH5, popt.printASCII);
                }


            }
            if (particle_index_end >= AllStreamlines.size() - 1){
                break;
            }
            particle_index_start = particle_index_end;
            iter++;
        }
        std::cout << "Total simulation Time : " << Simulation_time << std::endl;
    }

	void ParticleTrace::run_thread(int ithread) {
	    //if (ithread != 0){
	    //    bool stop= true;
	    //    std::cout << "Stop Here" << std::endl;
	    //}
        if (world.size() > 1 || (world.size() == 1 && popt.Nthreads == 1)){
            double show_every = popt.OutputFrequency;
            double cnt = 0.0;
            double Nsize = static_cast<double>(Streamlines4thread.size());
            for (int i = 0; i < Streamlines4thread.size(); ++i){
                if (popt.OutputFrequency > 0.1){
                    if (100*cnt/Nsize > show_every){
                        std::cout << "-" << show_every << " % |" << i << " of " << Nsize << std::endl;
                        show_every = show_every + popt.OutputFrequency;
                    }
                }
                //std::cout << ithread << ": " << i << std::endl;
                //Streamline S(Streamlines4thread[i].getEid(), Streamlines4thread[i].getSid(),Streamlines4thread[i].getLastParticle());
                bool tf;
                Domain.bisInProcessorPolygon(Streamlines4thread[i].getLastParticle().getP(),tf);
                if (!tf)
                    continue;

                ExitReason er = traceInner(Streamlines4thread[i]);
                Streamlines4thread[i].Close(er);
                cnt++;
            }
        }
        else{
            //TODO
            std::cout << "True multi threaded function is currently disabled" << std::endl;
        }

	}
}
