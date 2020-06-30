#pragma once
#include <iostream>
#include <fstream>
#include <string>

#include <boost/geometry/algorithms/assign.hpp>

#include "ichnos_structures.h"



namespace ICHNOS {

	namespace DEBUG {
		void displayParticleasVex(const Particle& P, bool prinAttr) {
			vec3 vn = P.getV().normalize();
			std::cout << std::setprecision(5) << std::fixed << "p = addpoint(0,{" << P.getP().x << "," << P.getP().z << "," << P.getP().y << "});";
			if (prinAttr)
				std::cout << std::setprecision(5) << std::fixed << " setpointattrib(0,'N',p,{" << vn.x << "," << vn.z << "," << vn.y << "},'set');";
			std::cout << std::endl;
		}
		void displayVectorasVex(vec3 p) {
			std::cout << std::setprecision(3) << std::fixed << "p = addpoint(0,{" << p.x << "," << p.z << "," << p.y << "});" << std::endl;
		}
		void displayPVasVex(vec3 p, vec3 v) {
			vec3 vn = v.normalize();
			std::cout << std::setprecision(3) << std::fixed << "p = addpoint(0,{" << p.x << "," << p.z << "," << p.y << "});";
			std::cout << std::setprecision(3) << std::fixed << " setpointattrib(0,'N',p,{" << vn.x << "," << vn.z << "," << vn.y << "},'set');";
			std::cout << std::endl;
		}
	}

	std::string num2Padstr(int i, int n) {
		std::stringstream ss;
		ss << std::setw(n) << std::setfill('0') << i;
		return ss.str();
	}

	void linspace(double min, double max, int n, std::vector<double>& v){
		v.clear();
		int iterator = 0;
		for (int i = 0; i <= n-2; i++){
			double temp = min + i*(max-min)/(floor(static_cast<double>(n)) - 1);
			v.insert(v.begin() + iterator, temp);
        	iterator += 1;
		}
		v.insert(v.begin() + iterator, max);
	}

	// https://stackoverflow.com/questions/2390912/checking-for-an-empty-file-in-c
	bool is_file_empty(std::ifstream& pFile){
		return pFile.peek() == std::ifstream::traits_type::eof();
	}

	void distributeParticlesAroundWellLayered(int eid, double x, double y, double top, double bot, std::vector<Streamline>& S, WellOptions wopt){
		std::vector<double> zval;
		linspace(bot, top, wopt.Nlayer, zval);
		int Nppl = wopt.Nparticles/wopt.Nlayer;
		double rads = (2.0*M_PI)/Nppl;
		std::vector<double> rads1; 
		linspace(0.0, 2.0*M_PI, wopt.Nlayer, rads1);
		std::vector<std::vector<double> > radpos;

		for (unsigned int i = 0; i < static_cast<unsigned int>(wopt.Nlayer); ++i){
			std::vector<double> tmp;
			linspace(0 + rads/2.0 + rads1[i],
								2.0*M_PI - rads/2.0 + rads1[i], Nppl, tmp);
			radpos.push_back(tmp);
		}
		
		int sid = 0;
		for(int i = 0; i < wopt.Nlayer; ++i){
			for( int j = 0; j < Nppl; ++j){
				vec3 temp;
				temp.x = cos(radpos[i][j])*wopt.Radius + x;
				temp.y = sin(radpos[i][j])*wopt.Radius + y;
				temp.z = zval[i];
				S.push_back(Streamline(eid, sid, Particle(temp)));
				//DEBUG::displayVectorasVex(temp);
				sid++;
			}
		}
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
		for (double i = 0; i <= 1+dt/2.0; i = i + dt){
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

		void readTopBot(std::string filename, PointSet2& Pset, bool top, bool bot) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
			}
			else {
				std::vector< std::pair<cgal_point_2,elev_data> > xy_data;
				elev_data el_data;
				std::string line;
				double x, y;
				while (getline(datafile, line)) {
					if (line.size() > 1) {
						std::istringstream inp(line.c_str());
						inp >> x;
						inp >> y;
						if (top && bot){
							inp >> el_data.top;
							inp >> el_data.bot;
						}
						else if (top){
							inp >> el_data.top;
						}
						else if (bot){
							inp >> el_data.bot;
						}
						cgal_point_2 p(x, y);
						xy_data.push_back(std::make_pair(p, el_data));
					}
				}
				Pset.insert(xy_data.begin(), xy_data.end());
				datafile.close();
			}
		}
		
		/*void read2DscalarField(std::string filename, pointCloud<double>& pntCld) {
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
		}*/

		void readNPSATVelocity(std::string filename, 
								std::vector<std::pair<cgal_point, NPSAT_data>>& data, 
								double VelocityMultiplier) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
			}
			else {
				std::string line;
				int ii;// , iproc;
				double x, y, z; // , vx, vy, vz;
				NPSAT_data npsat_data;
				ii = 0;
				while (getline(datafile, line)) {
					if (line.size() > 1) {
						std::istringstream inp(line.c_str());
						inp >> x;
						inp >> y;
						inp >> z;
						inp >> npsat_data.v.x;
						inp >> npsat_data.v.y;
						inp >> npsat_data.v.z;
						inp >> npsat_data.proc;
						npsat_data.v = npsat_data.v * VelocityMultiplier;
						npsat_data.id = ii;
						ii++;
						data.push_back(std::make_pair(cgal_point(x, y, z), npsat_data));
					}
				}
				std::cout << "Size of Data: " << data.size() << std::endl;
				datafile.close();
			}
		}
		
		void readNPSATVelocity(std::string filename, 
							std::vector<cgal_point_3>& pp,
							std::vector<NPSAT_data>& dd,
								double VelocityMultiplier) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
			}
			else {
				std::string line;
				int ii;// , iproc;
				double x, y, z; // , vx, vy, vz;
				NPSAT_data npsat_data;
				ii = 0;
				while (getline(datafile, line)) {
					if (line.size() > 1) {
						std::istringstream inp(line.c_str());
						inp >> x;
						inp >> y;
						inp >> z;
						inp >> npsat_data.v.x;
						inp >> npsat_data.v.y;
						inp >> npsat_data.v.z;
						inp >> npsat_data.proc;
						inp >> npsat_data.diameter;
						inp >> npsat_data.ratio;
						npsat_data.v = npsat_data.v * VelocityMultiplier;
						npsat_data.id = ii;
						ii++;
						//dd.push_back(std::make_pair(cgal_point_3(x, y, z), npsat_data));
						pp.push_back(cgal_point_3(x, y, z));
						dd.push_back(npsat_data);
					}
				}
				std::cout << "Size of Data: " << dd.size() << std::endl;
				datafile.close();
			}
		}
		

		/*void readVelocityFieldFile(std::string filename, pointCloud<vec3>& pntCld) {
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
						inp >> pntCld.NmaxPnts;
						inp >> pntCld.NminPnts;
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
		}*/

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
					distributeParticlesAroundWellLayered(eid, x, y, top, bot, S, wopt);
					//distributeParticlesAroundWell(eid, x, y, top, bot, S, wopt);
				}
				datafile.close();
				return true;
			}
		}

		bool readTrajfiles(std::string filename, streamlineMap& Smap) {
			std::ifstream datafile(filename.c_str());
			if (!datafile.good()) {
				std::cout << "Can't open the file" << filename << std::endl;
				return false;
			}
			else {
				if (is_file_empty(datafile)){
					return true;
				}
				std::string line, er_str;
				int eid, sid, pid, eid_prev, sid_prev;
				gPart particle;
				gStream streamline;
				streamlineMap::iterator eIt;
				std::map<int, gStream >::iterator sIt;
				eid_prev = -9;
				sid_prev = -9;
				while (getline(datafile, line)) {
					std::istringstream inp(line.c_str());
					inp >> pid;
					inp >> eid;
					inp >> sid;
					if (pid < 0) {
						inp >> er_str;
					}
					else {
						inp >> particle.p.x;
						inp >> particle.p.y;
						inp >> particle.p.z;
						inp >> particle.v.x;
						inp >> particle.v.y;
						inp >> particle.v.z;
						//inp >> particle.age;
					}

					bool dofind = true;
					if (eid_prev == eid && sid_prev == sid)
						dofind = false;

					if (dofind)
						eIt = Smap.find(eid);

					if (eIt != Smap.end()) {
						if (dofind)
							sIt = eIt->second.find(sid);
						if (sIt != eIt->second.end()) {
							if (pid < 0) {
								sIt->second.ex = castExitReasons2Enum(er_str);
							}
							else {
								sIt->second.particles.insert(std::pair<int, gPart>(pid, particle));
							}
						}
						else {
							gStream gs;
							if (pid < 0) {
								gs.ex = castExitReasons2Enum(er_str);
							}
							else {
								gs.particles.insert(std::pair<int, gPart>(pid, particle));
							}
							eIt->second.insert(std::pair<int, gStream >(sid, gs));
							sIt = eIt->second.find(sid);
						}
					}
					else {
						gStream gs;
						if (pid < 0) {
							gs.ex = castExitReasons2Enum(er_str);
						}
						else {
							gs.particles.insert(std::pair<int, gPart>(pid, particle));
						}
						std::map<int, gStream > tmpStream;
						tmpStream.insert(std::pair<int, gStream>(sid, gs));
						Smap.insert(std::pair<int, std::map<int, gStream > >(eid, tmpStream));
						// make sure that the iterators are pointing to the new data that were just inserted
						// most likely the following lines will belong to the same streamline
						eIt = Smap.find(eid);
						sIt = eIt->second.find(sid);
					}
					sid_prev = sid;
					eid_prev = eid;
				}
				return true;
			}
		}
	}

	

	namespace WRITE {
		void PrintParticle2Log(std::ofstream& log_file, Streamline& S, int i) {
			Particle pp = S.getParticle(i);
			log_file << pp.getPid() << " "
				<< S.getEid() << " "
				<< S.getSid() << " "
				<< std::setprecision(3) << std::fixed
				<< pp.getP().x << " " << pp.getP().y << " " << pp.getP().z << " "
				<< std::setprecision(4) << std::fixed
				<< pp.getV().x << " " << pp.getV().y << " " << pp.getV().z << " "
				/*<< pp.getAge()*/ << std::endl;
		}
		void PrintExitReason(std::ofstream& log_file, Streamline& S, ExitReason er) {
			log_file << -9 << " "
				<< S.getEid() << " "
				<< S.getSid() << " " 
				<< castExitReasons2String(er) << std::endl;
		}

		void printStreamslineMap(std::string filename, streamlineMap& SM) {
			std::ofstream out_file;
			out_file.open(filename.c_str());
			streamlineMap::iterator eit = SM.begin();
			std::map<int, gStream >::iterator sit;
			std::map<int, gPart>::iterator pit;
			
			for (; eit != SM.end(); ++eit) {
				sit = eit->second.begin();
				for (; sit != eit->second.end(); ++sit) {
					pit = sit->second.particles.begin();
					//double age = 0;
					//int i = 0;
					//vec3 p_prev, p_curr;
					//vec3 v_prev, v_curr, v_m;
					for (; pit != sit->second.particles.end(); ++pit) {
						/* // This is used to calculate the age
						// However to save some space in the output files
						// I wont print it
						if (i == 0) {
							p_prev = pit->second.p;
							v_prev = pit->second.v;
						}
						else {
							p_curr = pit->second.p;
							v_curr = pit->second.v;
							v_m = (v_prev + v_curr) * 0.5;
							double dst = (p_curr - p_prev).len();
							age += dst / v_m.len();
							p_prev = p_curr;
							v_prev = v_curr;
						}
						*/
						out_file << eit->first << " "
							<< sit->first << " "
							<< std::setprecision(2) << std::fixed
							<< pit->second.p.x << " "
							<< pit->second.p.y << " "
							<< pit->second.p.z << " "
							<< std::setprecision(5) << std::fixed
							<< pit->second.v.len() << " "
							// << pit->second.v.x << " "
							// << pit->second.v.y << " "
							// << pit->second.v.z << " " 
							/* << age */ << std::endl;
						//i++;
					}
					out_file << "-9 " 
							 << eit->first << " " 
							 << sit->first << " "
							 << castExitReasons2String(sit->second.ex)
							 << std::endl;
				}
			}
			out_file.close();
		}
	}

	void gather_particles(ICHNOS::options& opt) {
		for (int ireal = 0; ireal < opt.Popt.Nrealizations; ++ireal) {
			streamlineMap Smap;
			for (int i = 0; i < opt.niter; ++i) {
				for (int j = 0; j < opt.nproc; ++j) {
					std::string filename = opt.Popt.OutputFile + "_ireal_" + num2Padstr(ireal, 4) +
						"_iter_" + num2Padstr(i, 4) + "_proc_" + num2Padstr(j, 4) + ".traj";

					std::cout << "Reading ... " << filename << std::endl;
					bool tf = READ::readTrajfiles(filename, Smap);
					if (!tf) {
						std::cout << "Error while reading trajectory files" << std::endl;
						return;
					}
				}
				if (!opt.Popt.gatherOneFile) {
					std::string filename = opt.Popt.OutputFile + "_ireal_" + num2Padstr(ireal, 4) + "_gather_iter_" + num2Padstr(i, 4) + ".traj";
					ICHNOS::WRITE::printStreamslineMap(filename, Smap);
					Smap.clear();
				}
			}
			if (opt.Popt.gatherOneFile) {
				std::string filename = opt.Popt.OutputFile + "_ireal_" + num2Padstr(ireal, 4) + "_gather_ALL.traj";
				ICHNOS::WRITE::printStreamslineMap(filename, Smap);
			}
		}
	}

	//double interpolateScatter2D(cgal_Delaunay_2& T, coord_map& values, double x, double y) {
	//	cgal_point_2 p(x, y);
	//	std::vector<std::pair<cgal_point_2, coord_type>> coords;
	//	coord_type norm = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords)).second;
	//	coord_type res = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, value_access(values));
	//	return static_cast<double>(res);
	//}
	
	/*double interpolateScalarTree(std::unique_ptr<nano_kd_tree_scalar>& tree, vec3 p) {
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
					w = 1 / std::pow(std::sqrt(ret_matches[i].second), power);
					sumW += w;
					sumWVal += tree->dataset.kdtree_get_vel_vec(ret_matches[i].first) * w;
				}
			}
			return(sumWVal / sumW);
		}
	}*/

	//void interpolateIntTree(std::map<int, double>& id_dst, std::map<int, double>& proc_map,
	//	std::unique_ptr <nano_kd_tree_int>& tree, vec3 p) {
	//	double power = tree->dataset.Power;
	//	double threshold = tree->dataset.Threshold;
	//	double radius = tree->dataset.Radious;
	//	id_dst.clear();
	//	proc_map.clear();
	//	double query_pt[3] = { p.x, p.y, 0.0 };
	//	// This is a vector of pairs (id distance).  
	//	std::vector<std::pair<size_t, double> >   ret_matches;
	//	nanoflann::SearchParams params;
	//	size_t N = tree->radiusSearch(&query_pt[0], radius * radius, ret_matches, params);
	//	if (N == 0) {
	//		std::cerr << "There are no points around (" << p.x << "," << p.y << ") whithin " << radius << "distance" << std::endl;
	//		std::cerr << "consider increasing the radius" << std::endl;
	//	}
	//	else {
	//		double sumW = 0.0;
	//		double w;
	//		int iproc;
	//		std::map<int, double>::iterator it;
	//		for (size_t i = 0; i < N; i++) {
	//			iproc = tree->dataset.kdtree_get_proc(ret_matches[i].first);
	//			if (ret_matches[i].second < threshold) {
	//				id_dst.clear();
	//				id_dst.insert(std::pair<int, double>(static_cast<int>(ret_matches[i].first), 1.0));
	//				break;
	//			}
	//			else {
	//				w = 1 / std::pow(std::sqrt(ret_matches[i].second), power);
	//				it = proc_map.find(iproc);
	//				if (it == proc_map.end()) {
	//					proc_map.insert(std::pair<int, double>(iproc, w));
	//				}
	//				else {
	//					it->second += w;
	//				}
	//				sumW += w;
	//				id_dst.insert(std::pair<int, double>(static_cast<int>(ret_matches[i].first), w));
	//			}
	//		}
	//		it = id_dst.begin();
	//		for (it; it != id_dst.end(); ++it) {
	//			it->second = it->second / sumW;
	//		}
	//		it = proc_map.begin();
	//		for (it; it != proc_map.end(); ++it) {
	//			it->second = it->second / sumW;
	//		}
	//	}
	//}

	//void interpolateVectorTree(vec3& vel, std::map<int, double>& proc_map,
	//	std::unique_ptr <nano_kd_tree_vector>& tree, vec3 p, float scale = 1) {
	//	double power = tree->dataset.Power;
	//	double threshold = tree->dataset.Threshold;
	//	double radius = tree->dataset.Radious;
	//	size_t Nmax = static_cast<size_t>(tree->dataset.NmaxPnts);
	//	size_t Nmin = static_cast<size_t>(tree->dataset.NminPnts);
	//	proc_map.clear();
	//	double query_pt[3] = { p.x, p.y, p.z };

	//	std::vector<std::pair<size_t, double> >   ret_matches;
	//	nanoflann::SearchParams params;
	//	size_t N = 0;
	//	int count_attempts = 0;
	//	while (N < Nmin) {
	//		N = tree->radiusSearch(&query_pt[0], radius * radius, ret_matches, params);
	//		radius = radius * 2;
	//		count_attempts++;
	//		if (count_attempts > 20) {
	//			std::cout << "There are no points around (" << p.x << "," << p.y << "," << p.z << ") whithin " << radius << "  units radius" << std::endl;
	//			vel = vec3();
	//			return;
	//		}
	//	}
	//	N = std::min(N, Nmax);

	//	// find the largest horizontal and vertical distance between the point in question and the
	//	// closest points around this
	//	std::vector<vec3> closest_pnts(N);
	//	double xymax = 0;
	//	double zmax = 0;
	//	double tmpxy, tmpz;
	//	for (size_t i = 0; i < N; i++) {
	//		//std::cout << tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first).x << ", "
	//		//	<< tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first).y << ", "
	//		//	<< tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first).z << std::endl;
	//		closest_pnts[i] = tree->dataset.kdtree_get_pnt_vec(ret_matches[i].first) - p;
	//		tmpxy = std::sqrt(closest_pnts[i].x * closest_pnts[i].x + closest_pnts[i].y * closest_pnts[i].y);
	//		tmpz = std::sqrt(closest_pnts[i].z * closest_pnts[i].z);
	//		if (tmpxy > xymax)
	//			xymax = tmpxy;
	//		if (tmpz > zmax)
	//			zmax = tmpz;
	//	}

	//	double ratio = xymax / zmax;

	//	// interpolate using the scaled value
	//	double w, d;
	//	double sumW = 0;
	//	vec3 sumWVal;

	//	std::map<int, double>::iterator it;
	//	int iproc;

	//	for (size_t i = 0; i < N; i++) {
	//		closest_pnts[i].z = closest_pnts[i].z * ratio * scale + closest_pnts[i].z * (1.0 - scale);
	//		d = std::sqrt(closest_pnts[i].x * closest_pnts[i].x +
	//			closest_pnts[i].y * closest_pnts[i].y +
	//			(closest_pnts[i].z * closest_pnts[i].z));
	//		iproc = tree->dataset.kdtree_get_proc(ret_matches[i].first);
	//		if (d < threshold) {
	//			vel = tree->dataset.kdtree_get_vel_vec(ret_matches[i].first);
	//		}
	//		else {
	//			w = 1 / std::pow(d, power);
	//			it = proc_map.find(iproc);
	//			if (it == proc_map.end()) {
	//				proc_map.insert(std::pair<int, double>(iproc, w));
	//			}
	//			else {
	//				it->second += w;
	//			}
	//			sumW += w;
	//			sumWVal = sumWVal + tree->dataset.kdtree_get_vel_vec(ret_matches[i].first) * w;
	//		}
	//	}

	//	it = proc_map.begin();
	//	for (it; it != proc_map.end(); ++it) {
	//		it->second = it->second / sumW;
	//	}

	//	vel = sumWVal * (1 / sumW);
	//	
	//}

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
