#pragma once
#include <boost/mpi.hpp>
#include <random>

#include "ichnos_structures.h"
#include "velocity_base.h"
#include "ichnos_domain.h"


namespace ICHNOS {

	class ParticleTrace {
	public:
		ParticleTrace(
			boost::mpi::communicator& world_in,
			velocityField &VF_in,
			DomainBase &dom_in,
			ParticleOptions& popt_in);

		void Trace();

		bool readInputFiles(std::vector<Streamline>& S);

	private:
		const ParticleOptions popt;
		velocityField& VF;
		DomainBase& Domain;
		boost::mpi::communicator world;
		bool findNextPoint(const Particle& P, vec3& pnew, ExitReason& er);
		ExitReason traceInner(Streamline& S);
		void traceOuter(std::vector<Streamline>& S, int iter);
		ExitReason CheckNewPoint(vec3& p);
		ExitReason CheckNewPointAndCalcVelocity(vec3& p, vec3& v, int &proc);

		bool EulerStep(const Particle& P, vec3& pnew, ExitReason& er);
		bool RK2Step(const Particle& P, vec3& pnew, ExitReason& er);
		bool RK4Step(const Particle& P, vec3& pnew, ExitReason& er);
		bool RK45Step(const Particle& P, vec3& pnew, ExitReason& er);

		const std::vector<std::vector<double> > CF = getRK45coef();

		double adaptStepSize;
	};

	ParticleTrace::ParticleTrace(
		boost::mpi::communicator& world_in,
		velocityField& VF_in,
		DomainBase& dom_in,
		ParticleOptions &popt_in)
		:
		world(world_in),
		VF(VF_in),
		Domain(dom_in),
		popt(popt_in)
	{
		adaptStepSize = popt.StepSize;
	}

	bool ParticleTrace::readInputFiles(std::vector <Streamline>& S) {
		bool tfw = false;
		bool tfp = false;
		if (!popt.WellFile.empty()) {
			tfw = READ::readWellFile(popt.WellFile, S);
		}

		if (!popt.ParticleFile.empty()) {
			tfp = READ::readParticleFile(popt.ParticleFile, S);
		}

		if (tfw || tfp)
			return true;
		else
			return false;
	}

	void ParticleTrace::Trace() {
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

		int outer_iter = 0;
		while (true) {
			Part_Streamlines.clear();
			world.barrier();
			if (my_rank == 0) {
				std::default_random_engine generator;
				int N = popt.ParticlesInParallel;
				if (static_cast<int>(ALL_Streamlines.size()) < popt.ParticlesInParallel + 1000) {
					N = ALL_Streamlines.size();
				}

				for (int i = 0; i < N; ++i) {
					if (ALL_Streamlines.size() == 0)
						break;
					// It better to make this as singleton generator
					std::uniform_int_distribution<int> distribution(0, ALL_Streamlines.size() - 1);
					int ii = distribution(generator);
					Part_Streamlines.push_back(ALL_Streamlines[ii]);
					ALL_Streamlines.erase(ALL_Streamlines.begin() + ii);
				}
				//std::cout << "There are " << ALL_Streamlines.size() << " Streamlines to trace" << std::endl;
			}
			world.barrier();
			//std::cout << my_rank << " Before: " << Part_Streamlines.size() << std::endl;
			MPI::Sent_receive_Initial_streamlines_from0(Part_Streamlines, my_rank, world);
			//std::cout << my_rank << " After: " << Part_Streamlines.size() << std::endl;

			traceOuter(Part_Streamlines, outer_iter);
			
			break;
		}
	}

	void ParticleTrace::traceOuter(std::vector<Streamline>& S, int iter) {
		int my_rank = world.rank();
		int n_proc = world.size();

		const std::string log_file_name = (popt.OutputFile + "_iter_" + num2Padstr(iter, 4) + "_proc_" + num2Padstr(my_rank+1, 4) + ".traj");

		std::ofstream log_file;
		log_file.open(log_file_name.c_str());

		int trace_iter = 0;
		std::vector<Streamline> Snew;
		while (true) {
			Snew.clear();
			for (unsigned int i = 0; i < S.size(); ++i) {
				ExitReason er = traceInner(S[i]);

				if (er == ExitReason::FIRST_POINT_GHOST || er == ExitReason::FAR_AWAY)
					continue;

				for (unsigned int j = 0; j < S[i].size()-1; ++j) {
					WRITE::PrintParticle2Log(log_file, S[i], j);
				}

				if (er == ExitReason::CHANGE_PROCESSOR) {
					Snew.push_back(Streamline(
						S[i].getEid(), S[i].getSid(), S[i].getLastParticle(), S[i].getBBlow(), S[i].getBBupp(), S[i].StuckIter()));
				}
				else {
					WRITE::PrintParticle2Log(log_file, S[i], S[i].size() - 1);
					WRITE::PrintExitReason(log_file, S[i], er);
				}
			}

			world.barrier();
			std::cout << "Proc " << my_rank << " will send " << Snew.size() << " particles" << std::endl;
			int N_part2track = Snew.size();
			//std::cout << "Proc " << my_rank << " N_part2send = " << N_part2send << std::endl;
			MPI::sumScalar<int>(N_part2track, n_proc, world, MPI_INT);
			//std::cout << "Proc (after) " << my_rank << " N_part2send = " << N_part2send << std::endl;
			world.barrier();
			
			if (my_rank == 0)
				std::cout << "          Number of active particles: " << N_part2track << " --------" << std::endl << std::flush;

			if (N_part2track == 0)
				break;

			MPI::Send_receive_streamlines(Snew, S, world);
			std::cout << my_rank << " S new size: " << S.size() << std::endl;
			
			if (trace_iter == 5) {
				std::cout << "exit after " << trace_iter << " Iteration. DONT FORGET TO REMOVE THIS Condition" << std::endl;
				break;
			}

			trace_iter++;
			if (trace_iter > popt.MaxProcessorExchanges)
				break;
		}
		log_file.close();
	}

	ExitReason ParticleTrace::traceInner(Streamline& S) {
		vec3 v;
		vec3 p = S.getLastParticle().getP();
		int proc = world.rank();
		// Check if the particle is in the domain and exit imediately if not
		ExitReason er = CheckNewPoint(p);
		if (er == ExitReason::EXIT_SIDE || er == ExitReason::EXIT_BOTTOM || er == ExitReason::EXIT_TOP)
			return ExitReason::INIT_OUT;

		// Attempt to calculate the velocity
		std::map<int, double> proc_map;
		VF.calcVelocity(v, proc_map, p);
		if (v.isInvalid())
			return ExitReason::FAR_AWAY;
		proc = VF.calcProcID(proc_map);
		if (proc != world.rank())
			return ExitReason::FIRST_POINT_GHOST;

		// If the particle id is 0 then we have to assign velocity
		// If its not zero that means that the velocity has been calculated by the processor
		// that was tracking this streamline before end up in the currect processor
		if (S.getLastParticle().getPid() == 0)
			S.InitialVelocity(v);


		//S.getLastParticle().displayAsVEX(true);

		int count_iterations = 0;
		while (er == ExitReason::NO_EXIT) {
			bool foundPoint = findNextPoint(S.getLastParticle(), p, er);
			if (foundPoint) {
				er = CheckNewPointAndCalcVelocity(p, v, proc);
				S.AddParticle(Particle(p, v, proc, S.getLastParticle()));
				//S.getLastParticle().displayAsVEX(true);
				// In the unlike event that the velocity of the point is indeed zero
				// we can avoid unnessecary iterations
				if (v.isZero())
					return ExitReason::STUCK;
			}

			count_iterations++;
			if (count_iterations > popt.MaxIterationsPerStreamline)
				return ExitReason::MAX_INNER_ITER;
			if (S.StuckIter() > popt.StuckIterations)
				return ExitReason::STUCK;
		}
		S.Close(er);
		return er;
	}

	ExitReason ParticleTrace::CheckNewPointAndCalcVelocity(vec3& p, vec3& v, int& proc) {
		// This function checks only if the point is still inside the domain
		// It does not check if its in the velocity field owned by this processor
		ExitReason exitreason = CheckNewPoint(p);
		std::map<int, double> proc_map;
		// If the point is outside the domain we can still caluclate velocity if
		// the velocity field allows it. However the particle tracking will be terminated
		if (exitreason != ExitReason::NO_EXIT) {
			if (VF.InterpolateOutsideDomain) {
				VF.calcVelocity(v, proc_map, p);
			}
			else {
				v = vec3();
			}
			return exitreason;
		}
		else {
			VF.calcVelocity(v, proc_map, p);
			if (v.isInvalid())
				return ExitReason::FAR_AWAY;

			bool tf = VF.bIsInGhostArea(proc_map);
			if (tf) {
				proc = VF.calcProcID(proc_map);
				return ExitReason::CHANGE_PROCESSOR;
			}
			else {
				proc = world.rank();
				return ExitReason::NO_EXIT;
			}
		}
	}

	ExitReason ParticleTrace::CheckNewPoint(vec3& p) {
		bool bInDomain = true;
		Domain.bisInPolygon(p, bInDomain);
		if (!bInDomain)
			return ExitReason::EXIT_SIDE;

		double top, bottom;
		Domain.getTopBottomElevation(p, top, bottom);
		if (p.z < bottom)
			return ExitReason::EXIT_BOTTOM;
		if (p.z > top)
			return ExitReason::EXIT_TOP;

		return ExitReason::NO_EXIT;
	}

	bool ParticleTrace::findNextPoint(const Particle& P, vec3& pnew, ExitReason& er) {
		switch (popt.method)
		{
		case SolutionMethods::Euler:
			return EulerStep(P, pnew, er);
			break;
		case SolutionMethods::RK2:
			return RK2Step(P, pnew, er);
			break;
		case SolutionMethods::RK4:
			return RK4Step(P, pnew, er);
			break;
		case SolutionMethods::RK45:
			return RK45Step(P, pnew, er);
			break;
		default:
			return EulerStep(P, pnew, er);
			break;
		}
	}

	bool ParticleTrace::EulerStep(const Particle& P, vec3& pnew, ExitReason& er) {
		pnew = P.getV()* (1 / P.getV().len())* popt.StepSize* popt.Direction + P.getP();
		return true;
	}

	bool ParticleTrace::RK2Step(const Particle& P, vec3& pnew, ExitReason& er) {
		// We set this value to the current processor however this is not actually used
		// as it does not get out of this function
		// The actual processor will be calculated outside this function
		int proc = world.rank();
		// Find the new position using the velocity of the current point
		vec3 p1 = P.getV()* (1 / P.getV().len())*popt.StepSize * popt.Direction + P.getP();
		// Find the velocity on that point
		vec3 v1;
		er = CheckNewPointAndCalcVelocity(p1, v1, proc);
		if (v1.isZero()) return false;
		// Take a step with the average velocity
		vec3 v_m = v1 * 0.5 + P.getV() * 0.5;
		pnew = P.getV()* (1 / v_m.len())* popt.StepSize* popt.Direction + P.getP();
		return true;
	}

	bool ParticleTrace::RK4Step(const Particle& P, vec3& pnew, ExitReason& er) {
		// See RK2Step
		int proc = world.rank();
		// Step 1 ----------------------------
		// take a half step from the current point using the velocity of that point
		vec3 p2 = P.getV() * (1 / (P.getV().len()/2)) * popt.StepSize * popt.Direction + P.getP();
		// Find the velocity on point p2 that point
		vec3 v2;
		er = CheckNewPointAndCalcVelocity(p2, v2, proc);
		if (v2.isZero()) return false;

		// Step 2 ----------------------------
		// take a half step from the current point using the v2 velocity
		vec3 p3 = v2 * (1 / (v2.len() / 2)) * popt.StepSize * popt.Direction + P.getP();
		// Find the velocity on point p3 that point
		vec3 v3;
		er = CheckNewPointAndCalcVelocity(p3, v3, proc);
		if (v3.isZero()) return false;

		// Step 3 ----------------------------
		// take a full step from the current point using the v3 velocity
		vec3 p4 = v3 * (1 / v3.len()) * popt.StepSize * popt.Direction + P.getP();
		// Find the velocity on point p4 that point
		vec3 v4;
		er = CheckNewPointAndCalcVelocity(p4, v4, proc);
		if (v4.isZero()) return false;

		vec3 v_m = (P.getV() + v2 * 2 + v3 * 2 + v4)*(1/6);
		pnew = P.getV() * (1 / v_m.len()) * popt.StepSize * popt.Direction + P.getP();
		return true;
	}

	bool ParticleTrace::RK45Step(const Particle& P, vec3& pnew, ExitReason& er) {
		// See RK2Step
		int proc = world.rank();
		vec3 v_m, p2, p3, p4, p5, p6, v2, v3, v4, v5, v6, yn, zn;
		vec3 v1 = P.getV();

		double ssdir = adaptStepSize * popt.Direction;

		//Runge-Kutta-Fehlberg Method (RKF45)		// Step 1 ----------------------------		int i = 0;		p2 = v1 * ((1 / (v1.len()*CF[i][0])) * ssdir) + P.getP();
		er = CheckNewPointAndCalcVelocity(p2, v2, proc);
		if (v2.isZero()) { return false; }		// Step 2 ----------------------------		++i;		v_m = v1 * CF[i][1] + v2 * CF[i][2];		p3 = v_m*(1/v_m.len()*CF[i][0]) * ssdir + P.getP();
		er = CheckNewPointAndCalcVelocity(p3, v3, proc);
		if (v3.isZero()) return false;		// Step 3 ----------------------------		++i;		v_m = v1 * CF[i][1] + v2 * CF[i][2] + v3 * CF[i][3];		p4 = v_m * (1 / v_m.len() * CF[i][0]) * ssdir + P.getP();
		er = CheckNewPointAndCalcVelocity(p4, v4, proc);
		if (v4.isZero()) return false;		// Step 4 ----------------------------		++i;		v_m = v1 * CF[i][1] + v2 * CF[i][2] + v3 * CF[i][3] + v4 * CF[i][4];		p5 = v_m * (1 / v_m.len() * CF[i][0]) * ssdir + P.getP();
		er = CheckNewPointAndCalcVelocity(p5, v5, proc);
		if (v5.isZero()) return false;		// Step 5 ----------------------------		++i;		v_m = v1 * CF[i][1] + v2 * CF[i][2] + v3 * CF[i][3] + v4 * CF[i][4] + v5 * CF[i][5];		p6 = v_m * (1 / v_m.len() * CF[i][0]) * ssdir + P.getP();
		er = CheckNewPointAndCalcVelocity(p6, v6, proc);
		if (v6.isZero()) return false;		//// Calculate the point using two different paths		++i;		v_m = v1 * CF[i][1] + v3 * CF[i][2] + v4 * CF[i][3] + v5 * CF[i][4];		yn = v_m * (1 / v_m.len() * CF[i][0]) * ssdir + P.getP();		i++;		v_m = v1 * CF[i][1] + v3 * CF[i][2] + v4 * CF[i][3] + v5 * CF[i][4] + v6 * CF[i][5];		zn = v_m * (1 / v_m.len() * CF[i][0]) * ssdir + P.getP();		double R = (yn - zn).len();		if (R < popt.ToleranceStepSize) {			// We can increase the step size			adaptStepSize = 0.84 * pow(popt.ToleranceStepSize / R, 0.25) * popt.StepSize;			if (adaptStepSize > popt.MaxStepSize)				adaptStepSize = popt.MaxStepSize;			if (adaptStepSize < popt.MinStepSize)				adaptStepSize = popt.MinStepSize;			pnew = yn;		}		else {			if (adaptStepSize > popt.MinStepSize) {				adaptStepSize = 0.84 * pow(popt.ToleranceStepSize / R, 0.25) * popt.StepSize;				RK45Step(P, pnew, er);			}			else {				pnew = yn;				std::cout << "The algorithm will continue although it cannot reach the desired accuracy of "					<< popt.ToleranceStepSize					<< "with the limit of minimum step size " << popt.MinStepSize << std::endl;				std::cout << "Either relax the tolerance or reduce the Minimum Step Size" << std::endl;			}		}
		return true;
	}

}
