#pragma once
#include <iostream>
#include <fstream>
#include <string>

#include <boost/geometry/algorithms/assign.hpp>

#include "ichnos_structures.h"



namespace ICHNOS {

	namespace DEBUG {
		void displayParitcleasVex(const Particle& P, bool prinAttr) {
			vec3 vn = P.getV() * (1 / P.getV().len());
			std::cout << "p = addpoint(0,{" << P.getP().x << "," << P.getP().z << "," << P.getP().y << "});";
			if (prinAttr)
				std::cout << " setpointattrib(0,'N',p,{" << vn.x << "," << vn.z << "," << vn.y << "},'set');";
			std::cout << std::endl;
		}
		void displayVectorasVex(vec3 p) {
			std::cout << "p = addpoint(0,{" << p.x << "," << p.z << "," << p.y << "});" << std::endl;
		}
		void displayPVasVex(vec3 p, vec3 v) {
			vec3 vn = v * (1 / v.len());
			std::cout << "p = addpoint(0,{" << p.x << "," << p.z << "," << p.y << "});";
			std::cout << " setpointattrib(0,'N',p,{" << vn.x << "," << vn.z << "," << vn.y << "},'set');";
			std::cout << std::endl;
		}
	}

	std::string num2Padstr(int i, int n) {
		std::stringstream ss;
		ss << std::setw(n) << std::setfill('0') << i;
		return ss.str();
	}
	
	void distributeParticlesAroundWell(int eid, double x, double y, double top, double bot, std::vector<Streamline>& S, WellOptions wopt) {
		const double PI = 4 * atan(1);
		double maxt = 2*PI*static_cast<double>(wopt.Nlayer);
		double dt = 1 / (static_cast<double>(wopt.Nparticles) - 1);
		double t;
		double r = wopt.Radius;
		double h = top - bot;
		vec3 pp;

		int sid = 0;
		for (double i = 0; i <= 1; i = i + dt){
			t = i * maxt;
			pp.x = x + r * cos(t);
			pp.y = y + r * sin(t);
			pp.z = bot + h*i;
			//Particle(pp).displayAsVEX(false);
			S.push_back(Streamline(eid, sid, Particle(pp)));
			sid++;
		}
	}

	namespace READ {
		
		void read2DscalarField(std::string filename, pointCloud<double>& pntCld) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
			}
			else {
				std::string line;
				double x, y, z;
				getline(datafile, line);
				{
					std::istringstream inp(line.c_str());
					inp >> pntCld.Radious;
					inp >> pntCld.Power;
					inp >> pntCld.Threshold;
				}
				while (getline(datafile, line)) {
					std::istringstream inp(line.c_str());
					inp >> x;
					inp >> y;
					inp >> z;
					pntCld.pts.push_back(vec3(x, y, 0.0));
					pntCld.vel.push_back(z);
				}
				datafile.close();
			}
		}

		void readVelocityFieldFile(std::string filename, pointCloud<vec3>& pntCld) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
			}
			else {
				std::string comment("#");
				std::string line;
				int ii;
				double x, y, z;
				int readSection = 0;
				while (getline(datafile, line)) {
					if (line.size() > 1)
						if (comment.compare(0, 1, line, 0, 1) == 0)
							continue;
					if (readSection == 0) {
						std::istringstream inp(line.c_str());
						inp >> pntCld.Radious;
						inp >> pntCld.Power;
						inp >> pntCld.Threshold;
						inp >> ii;
						if (ii == 0)
							pntCld.InterpolateOutside = false;
						else
							pntCld.InterpolateOutside = true;
						readSection++;
						continue;
					}

					if (readSection == 1) {
						std::istringstream inp(line.c_str());
						inp >> ii;
						inp >> x;
						inp >> y;
						inp >> z;
						pntCld.proc.push_back(ii);
						pntCld.pts.push_back(vec3(x, y, z));
						inp >> x;
						inp >> y;
						inp >> z;
						pntCld.vel.push_back(vec3(x, y, z));
					}
				}
				datafile.close();

			}
		}

		void readPolygonDomain(std::string filename, boostPolygon& polygon) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
			}
			else {
				std::vector<boostPoint> polygonPoints;
				std::string line;
				double x, y;
				while (getline(datafile, line)) {
					std::istringstream inp(line.c_str());
					inp >> x;
					inp >> y;
					polygonPoints.push_back(boostPoint(x, y));
				}
				datafile.close();

				boost::geometry::assign_points(polygon, polygonPoints);
				boost::geometry::correct(polygon);
			}
		}

		bool readParticleFile(std::string filename, std::vector <Streamline>& S) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
				return false;
			}
			else {
				std::string comment("#");
				std::string line;
				int eid, sid;
				vec3 pp;
				while (getline(datafile, line)) {
					if (line.size() > 1)
						if (comment.compare(0, 1, line, 0, 1) == 0)
							continue;

					std::istringstream inp(line.c_str());
					inp >> eid;
					inp >> sid;
					inp >> pp.x;
					inp >> pp.y;
					inp >> pp.z;
					//P.push_back(Particle<T>(pp));
					S.push_back(Streamline(eid, sid, Particle(pp)));
				}
				datafile.close();
				return true;
			}
		}

		bool readWellFile(std::string filename, std::vector <Streamline>& S) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
				return false;
			}
			else {
				WellOptions wopt;
				std::string line;
				int eid;
				double x, y, top, bot;
				getline(datafile, line);
				{
					std::istringstream inp(line.c_str());
					inp >> wopt.Nparticles;
					inp >> wopt.Nlayer;
					inp >> wopt.Radius;
				}
				while (getline(datafile, line)) {
					std::istringstream inp(line.c_str());
					inp >> eid;
					inp >> x;
					inp >> y;
					inp >> top;
					inp >> bot;
					distributeParticlesAroundWell(eid, x, y, top, bot, S, wopt);
				}
				datafile.close();
				return true;
			}
		}
	}

	namespace WRITE {
		void PrintParticle2Log(std::ofstream& log_file, Streamline& S, int i) {
			Particle pp = S.getParticle(i);
			log_file << pp.getPid() << " \t"
				<< S.getEid() << " \t"
				<< S.getSid() << " \t"
				<< std::setprecision(3) << std::fixed
				<< pp.getP().x << " \t" << pp.getP().y << " \t" << pp.getP().z << " \t"
				<< std::setprecision(4) << std::fixed
				<< pp.getV().x << " \t" << pp.getV().y << " \t" << pp.getV().z << " \t"
				<< pp.getAge() << std::endl;
		}
		void PrintExitReason(std::ofstream& log_file, Streamline& S, ExitReason er) {
			log_file << -9 << " \t"
				<< S.getEid() << " \t"
				<< S.getSid() << " \t" 
				<< castExitReasons2String(er) << std::endl;
		}
	}
	
	double interpolateScalarTree(std::unique_ptr<nano_kd_tree_scalar>& tree, vec3 p) {
		double thres = tree->dataset.Threshold;
		double radius = tree->dataset.Radious;
		double power = tree->dataset.Power;
		double query_pt[3] = { p.x, p.y, 0.0 };
		std::vector<std::pair<size_t, double> >   ret_matches;
		nanoflann::SearchParams params;
		size_t N = tree->radiusSearch(&query_pt[0], radius*radius, ret_matches, params);

		if (N == 0) {
			std::cerr << "There are no points around (" << p.x << "," << p.y << ") whithin " << radius << "distance" << std::endl;
			std::cerr << "consider increasing the radius" << std::endl;
			return -999999;
		}
		else {
			double sumW = 0;
			double sumWVal = 0;
			double w;
			for (size_t i = 0; i < N; i++) {
				if (ret_matches[i].second < thres) {
					return tree->dataset.kdtree_get_pt(ret_matches[i].first, 3);
				}
				else {
					w = 1 / std::pow(ret_matches[i].second, power);
					sumW += w;
					sumWVal += tree->dataset.kdtree_get_pt(ret_matches[i].first, 3) * w;
				}
			}
			return(sumWVal / sumW);
		}
	}

	void interpolateVectorTree(vec3& vel, std::map<int, double>& proc_map,
		std::unique_ptr <nano_kd_tree_vector>& tree, vec3 p, float scale = 1) {
		double power = tree->dataset.Power;
		double threshold = tree->dataset.Threshold;
		double radius = tree->dataset.Radious;
		proc_map.clear();
		double query_pt[3] = { p.x, p.y, p.z };

		std::vector<std::pair<size_t, double> >   ret_matches;
		nanoflann::SearchParams params;
		size_t N = tree->radiusSearch(&query_pt[0], radius * radius, ret_matches, params);
		if (N == 0) {
			//std::cerr << "There are no points around (" << p.x << "," << p.y << "," << p.z << ") whithin " << radius << "  units radius" << std::endl;
			//std::cerr << "consider increasing the radius" << std::endl;
			vel = vec3(-99999, -99999, -99999);
		}
		else {
			// find the largest horizontal and vertical distance between the point in question and the
			// closest points around this
			std::vector<vec3> closest_pnts(N);
			double xymax = 0;
			double zmax = 0;
			double tmpxy, tmpz;

			for (size_t i = 0; i < N; i++) {
				//std::cout << tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first).x << ", "
				//	<< tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first).y << ", "
				//	<< tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first).z << std::endl;
				closest_pnts[i] = tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first) - p;
				tmpxy = std::sqrt(closest_pnts[i].x * closest_pnts[i].x + closest_pnts[i].y * closest_pnts[i].y);
				tmpz = std::sqrt(closest_pnts[i].z * closest_pnts[i].z);
				if (tmpxy > xymax)
					xymax = tmpxy;
				if (tmpz > zmax)
					zmax = tmpz;
			}

			double ratio = xymax / zmax;

			// interpolate using the scaled value
			double w, d;
			double sumW = 0;
			vec3 sumWVal;

			std::map<int, double>::iterator it;
			int iproc;

			for (size_t i = 0; i < N; i++) {
				closest_pnts[i].z = closest_pnts[i].z * ratio * scale + closest_pnts[i].z * (1 - scale);
				d = std::sqrt(closest_pnts[i].x * closest_pnts[i].x +
					closest_pnts[i].y * closest_pnts[i].y +
					(closest_pnts[i].z * closest_pnts[i].z));
				iproc = tree->dataset.kdtree_get_proc(ret_matches[i].first);
				if (d < threshold) {
					vel = tree->dataset.kdtree_get_vel_vec(ret_matches[i].first);
				}
				else {
					w = 1 / std::pow(d, power);
					it = proc_map.find(iproc);
					if (it == proc_map.end()) {
						proc_map.insert(std::pair<int, double>(iproc, w));
					}
					else {
						it->second += w;
					}
					sumW += w;
					sumWVal = sumWVal + tree->dataset.kdtree_get_vel_vec(ret_matches[i].first) * w;
				}
			}

			it = proc_map.begin();
			for (it; it != proc_map.end(); ++it) {
				it->second = it->second / sumW;
			}

			vel = sumWVal * (1 / sumW);
		}
	}

	bool is_input_scalar(std::string input) {
		// try to convert the input to scalar
		bool outcome;
		try {
			double value = std::stod(input);
			value++; // something to surpress the warning
			outcome = true;
		}
		catch (...) {
			outcome = false;
		}
		return outcome;
	}

	namespace MPI {
		namespace DEBUG {
			template <typename T>
			void DataPerProc(std::vector<T>& v, int myrank, std::string vname) {
				for (unsigned int i = 0; i < v.size(); ++i)
					std::cout << "Proc: " << myrank << ": " << vname << "[" << i << "] = " << v[i] << std::endl;
			}

			template <typename T>
			void PrintVectorData(std::vector<std::vector<T> >& v, int myrank, std::string vname) {
				for (int i = 0; i < v.size(); i++){
					for (int j = 0; j < v[i].size(); ++j) {
						std::cout << "Proc: " << myrank << ": " << vname << "[" << i << "][" << j << "] = " << v[i][j] << std::endl;
					}
				}
			}

			template<typename T>
			void PrintVectorSize(std::vector<std::vector<T> >& v, int myrank, std::string vname) {
				if (myrank < static_cast<int>(v.size()))
					std::cout << "I'm proc " << myrank << " and I have " << v[myrank].size() << " elements of " << vname << std::endl;
				else
					std::cout << "Vector " << vname << " doesnt have space for rank " << myrank << std::endl;

			}
		}
		

		

		void Send_receive_size(int N, int n_proc, std::vector<int>& output, boost::mpi::communicator& world) {

			output.clear();
			output.resize(n_proc);
			std::vector<int> temp(n_proc, 1);
			std::vector<int> displs(n_proc);
			for (int i = 1; i < n_proc; i++)
				displs[i] = displs[i - 1] + 1;

			MPI_Allgatherv(&N, // This is what this processor will send to every other
				1, //This is the size of the message from this processor
				MPI_INT, // The data type will be sent
				&output[0], // This is where the data will be send on each processor
				&temp[0], // an array with the number of points to be sent/receive
				&displs[0],
				MPI_INT, // The data type will be received
				world);

			//for (int i = 0; i < n_proc; ++i)
			//	std::cout << world.rank() << ":" << output[i] << std::endl;
		}

		template <typename T>
		void Send_receive_data(std::vector<std::vector<T> >& data,
			std::vector <int> N_data_per_proc,
			unsigned int my_rank,
			boost::mpi::communicator& world,
			MPI_Datatype MPI_TYPE) {

			// data is a vector of vectors of type T with size equal to n_proc.
			// This function transfer to all processors the content of data[my_rank]
			// if there are any data in data[i], where i=[1,n_proc; i!=myrank] this will be deleted
			// The size of data[my_rank].size() = N_data_per_proc[my_rank]. This is user responsibility

			int N = data[my_rank].size();
			if (N == 0)
				data[my_rank].push_back(0);

			unsigned int n_proc = data.size();
			std::vector<int> displs(n_proc);
			displs[0] = 0;
			for (unsigned int i = 1; i < n_proc; i++)
				displs[i] = displs[i - 1] + N_data_per_proc[i - 1];
			int totdata = displs[n_proc - 1] + N_data_per_proc[n_proc - 1];
			
			std::vector<T> temp_receive(totdata);

			
			MPI_Allgatherv(&data[my_rank][0], // This is what this processor will send to every other
				N, //This is the size of the message from this processor
				MPI_TYPE, // The data type will be sent
				&temp_receive[0], // This is where the data will be send on each processor
				&N_data_per_proc[0], // an array with the number of points to be sent/receive
				&displs[0],
				MPI_TYPE, world);

			//for (int i = 0; i < temp_receive.size(); ++i)
			//	std::cout << world.rank() << ": temp_receive[" << i << "] = " << temp_receive[i] << std::endl;

			
			// Now put the data in the data vector
			for (unsigned int i = 0; i < n_proc; ++i) {
				data[i].clear();
				data[i].resize(N_data_per_proc[i]);
				for (int j = 0; j < N_data_per_proc[i]; ++j)
					data[i][j] = temp_receive[displs[i] + j];
			}
		}

		template <typename T>
		void broadcast_vector (std::vector<T> & v, int proc, boost::mpi::communicator& world){
			T value;
			int VectorSize;

			if (world.rank() == proc) 
				VectorSize = static_cast<int>(v.size());

			boost::mpi::broadcast(world, VectorSize, proc);
			if (world.rank() != proc)
				v.clear();

			for (int i = 0; i < VectorSize; ++i) {
				if (world.rank() == proc) {
					value = v[i];
				}
				boost::mpi::broadcast(world, value, proc);
				if (world.rank() != proc) {
					v.push_back(value);
				}
			}
		}

		template <typename T>
		void sumScalar(T& scalar, int n_proc, boost::mpi::communicator& world, MPI_Datatype MPI_TYPE) {
			std::vector<T> Allscalar(n_proc);
			std::vector<int> temp(n_proc, 1);
			std::vector<int> displs(n_proc);
			for (int i = 1; i < n_proc; i++)
				displs[i] = displs[i - 1] + 1;

			MPI_Allgatherv(&scalar, // This is what this processor will send to every other
				1, //This is the size of the message from this processor
				MPI_TYPE, // The data type will be sent
				&Allscalar[0], // This is where the data will be send on each processor
				&temp[0], // an array with the number of points to be sent/receive
				&displs[0],
				MPI_TYPE, // The data type will be received
				world);
			scalar = 0;
			for (int i = 0; i < n_proc; ++i) {
				scalar += Allscalar[i];
			}
		}

		template <typename T>
		void maxScalar(T& scalar, int n_proc, boost::mpi::communicator& world, MPI_Datatype MPI_TYPE) {
			std::vector<T> Allscalar(n_proc);
			std::vector<int> temp(n_proc, 1);
			std::vector<int> displs(n_proc);
			for (unsigned int i = 1; i < n_proc; i++)
				displs[i] = displs[i - 1] + 1;

			MPI_Allgatherv(&scalar, // This is what this processor will send to every other
				1, //This is the size of the message from this processor
				MPI_TYPE, // The data type will be sent
				&Allscalar[0], // This is where the data will be send on each processor
				&temp[0], // an array with the number of points to be sent/receive
				&displs[0],
				MPI_TYPE, // The data type will be received
				world);
			scalar = -9999999999;
			for (int i = 0; i < n_proc; ++i) {
				if (Allscalar[i] > scalar)
					scalar = Allscalar[i];
			}
		}

		void Sent_receive_Initial_streamlines_from0(
			std::vector<Streamline>& S,
			int my_rank, boost::mpi::communicator& world) {
			world.barrier();
			int n_proc = world.size();
			if (n_proc == 1)
				return;

			std::vector<std::vector<double> > px(n_proc);
			std::vector<std::vector<double> > py(n_proc);
			std::vector<std::vector<double> > pz(n_proc);
			std::vector<std::vector<int> > E_id(n_proc);
			std::vector<std::vector<int> > S_id(n_proc);
			std::vector<std::vector<int> > proc_id(n_proc);

			// copy the data
			for (unsigned int i = 0; i < S.size(); ++i) {
				px[my_rank].push_back(S[i].getLastParticle().getP().x);
				py[my_rank].push_back(S[i].getLastParticle().getP().y);
				pz[my_rank].push_back(S[i].getLastParticle().getP().z);
				E_id[my_rank].push_back(S[i].getEid());
				S_id[my_rank].push_back(S[i].getSid());
			}
			world.barrier();
			// Send everything to every processor
			std::vector<int> data_per_proc;
			
			Send_receive_size(static_cast<int>(px[my_rank].size()), n_proc, data_per_proc, world);
			//DEBUG::DataPerProc<int>(data_per_proc, my_rank, "After");
			Send_receive_data<double>(px, data_per_proc, my_rank, world, MPI_DOUBLE);
			Send_receive_data<double>(py, data_per_proc, my_rank, world, MPI_DOUBLE);
			Send_receive_data<double>(pz, data_per_proc, my_rank, world, MPI_DOUBLE);

			Send_receive_data<int>(E_id, data_per_proc, my_rank, world, MPI_INT);
			Send_receive_data<int>(S_id, data_per_proc, my_rank, world, MPI_INT);
			//DEBUG::PrintVectorData<int>(S_id, my_rank, "E_id");
			S.clear();

			for (unsigned int i = 0; i < px.size(); ++i) 
				for (unsigned int j = 0; j < px[i].size(); ++j) 
					S.push_back(Streamline(E_id[i][j], S_id[i][j], Particle(vec3(px[i][j], py[i][j], pz[i][j]))));
			world.barrier();
		 }

		void Send_receive_streamlines(std::vector<Streamline>& Ssend, 
			std::vector<Streamline>& Srecv, boost::mpi::communicator& world) {
			int my_rank = world.rank();
			int n_proc = world.size();
			Srecv.clear();

			std::vector<std::vector<double> > px(n_proc);
			std::vector<std::vector<double> > py(n_proc);
			std::vector<std::vector<double> > pz(n_proc);
			std::vector<std::vector<double> > vx(n_proc);
			std::vector<std::vector<double> > vy(n_proc);
			std::vector<std::vector<double> > vz(n_proc);
			std::vector<std::vector<int> > E_id(n_proc);
			std::vector<std::vector<int> > S_id(n_proc);
			std::vector<std::vector<int> > proc_id(n_proc);
			std::vector<std::vector<int> > p_id(n_proc);
			std::vector<std::vector<int> > Nstuck(n_proc);
			std::vector<std::vector<double> > BBlx(n_proc);
			std::vector<std::vector<double> > BBly(n_proc);
			std::vector<std::vector<double> > BBlz(n_proc);
			std::vector<std::vector<double> > BBux(n_proc);
			std::vector<std::vector<double> > BBuy(n_proc);
			std::vector<std::vector<double> > BBuz(n_proc);
			world.barrier();
			for (unsigned int i = 0; i < Ssend.size(); ++i) {
				px[my_rank].push_back(Ssend[i].getLastParticle().getP().x);
				py[my_rank].push_back(Ssend[i].getLastParticle().getP().y);
				pz[my_rank].push_back(Ssend[i].getLastParticle().getP().z);
				vx[my_rank].push_back(Ssend[i].getLastParticle().getV().x);
				vy[my_rank].push_back(Ssend[i].getLastParticle().getV().y);
				vz[my_rank].push_back(Ssend[i].getLastParticle().getV().z);
				E_id[my_rank].push_back(Ssend[i].getEid());
				S_id[my_rank].push_back(Ssend[i].getSid());
				proc_id[my_rank].push_back(Ssend[i].getLastParticle().getProc());
				p_id[my_rank].push_back(Ssend[i].getLastParticle().getPid());
				Nstuck[my_rank].push_back(Ssend[i].StuckIter());
				BBlx[my_rank].push_back(Ssend[i].getBBlow().x);
				BBly[my_rank].push_back(Ssend[i].getBBlow().y);
				BBlz[my_rank].push_back(Ssend[i].getBBlow().z);
				BBux[my_rank].push_back(Ssend[i].getBBupp().x);
				BBuy[my_rank].push_back(Ssend[i].getBBupp().y);
				BBuz[my_rank].push_back(Ssend[i].getBBupp().z);
				
				
				//if (my_rank == 1) {
				//	std::cout << E_id[my_rank][i] << ", " << S_id[my_rank][i] << ", " << p_id[my_rank][i] << ", " << proc_id[my_rank][i] << ", "
				//		<< px[my_rank][i] << ", " << py[my_rank][i] << ", " << pz[my_rank][i] << ", "
				//		<< vx[my_rank][i] << ", " << vy[my_rank][i] << ", " << vz[my_rank][i] << std::endl;
				//}
				
			}
			world.barrier();
			std::vector<int> data_per_proc;
			MPI::Send_receive_size(px[my_rank].size(), n_proc, data_per_proc, world);
			//if (my_rank == 1) {
			//	for (int i = 0; i < data_per_proc.size(); ++i)
			//		std::cout << "data_per_proc " << data_per_proc[i] << std::endl;
			//}

			MPI::Send_receive_data<double>(px, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(py, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(pz, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(vx, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(vy, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(vz, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<int>(E_id, data_per_proc, my_rank, world, MPI_INT);
			MPI::Send_receive_data<int>(S_id, data_per_proc, my_rank, world, MPI_INT);
			MPI::Send_receive_data<int>(proc_id, data_per_proc, my_rank, world, MPI_INT);
			MPI::Send_receive_data<int>(p_id, data_per_proc, my_rank, world, MPI_INT);
			MPI::Send_receive_data<int>(Nstuck, data_per_proc, my_rank, world, MPI_INT);
			MPI::Send_receive_data<double>(BBlx, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(BBly, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(BBlz, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(BBux, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(BBuy, data_per_proc, my_rank, world, MPI_DOUBLE);
			MPI::Send_receive_data<double>(BBuz, data_per_proc, my_rank, world, MPI_DOUBLE);

			/*
			if (my_rank == 0) {
				for (int i = 0; i < px.size(); ++i) {
					for (int j = 0; j < px[i].size(); ++j) {
						std::cout << E_id[i][j] << ":" << S_id[i][j] << ":" << p_id[i][j] << ":" << proc_id[i][j] << std::endl;
						std::cout << "px[" << i << "][" << j << "]=" << px[i][j] <<  std::endl;
						std::cout << "py[" << i << "][" << j << "]=" << py[i][j] << std::endl;
						std::cout << "pz[" << i << "][" << j << "]=" << pz[i][j] << std::endl;

					}
				}
			}
			*/
			
			// Now loop through the data pick the ones that the other processors found 
			// that they should go on the current processor
			world.barrier();
			for (unsigned int i = 0; i < px.size(); ++i) {
				if (i == my_rank)
					continue;
				for (unsigned int j = 0; j < px[i].size(); ++j) {
					if (proc_id[i][j] == my_rank) {
						vec3 p = vec3(px[i][j], py[i][j], pz[i][j]);
						vec3 v = vec3(vx[i][j], vy[i][j], vz[i][j]);
						vec3 bl = vec3(BBlx[i][j], BBly[i][j], BBlz[i][j]);
						vec3 bu = vec3(BBux[i][j], BBuy[i][j], BBuz[i][j]);
						Streamline Stmp = Streamline(E_id[i][j], S_id[i][j],
							Particle(p, v, p_id[i][j], proc_id[i][j]),
							bl, bu, Nstuck[i][j]);
						Srecv.push_back(Stmp);
					}
				}
			}
			world.barrier();
			//std::cout << my_rank << " Srecv size: " << Srecv.size() << std::endl;
		}
	 } 
}
