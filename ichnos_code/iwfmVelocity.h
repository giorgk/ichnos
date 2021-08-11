#pragma once

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>

#include "ichnos_structures.h"
#include "velocity_base.h"

namespace po = boost::program_options;
namespace ic = ICHNOS;

namespace IWFM {
	class mat2 {
	public:
		mat2() {}
		void init(int nrow, int ncol){ M.resize(nrow, std::vector<double>(ncol)); }
		double& operator()(int i, int j) {
			return M[i][j];
		}
		void set(int i, int j, double v) {
			M[i][j] = v;
		}
	private:
		std::vector<std::vector<double> > M;
	};

	class xynode {
	public:
		xynode(){}
		xynode(ic::vec3 p)
			:
			pos(p)
		{}

		void init(int nrow, int ncol) {
			vx.init(nrow, ncol);
			vy.init(nrow, ncol);
			vz.init(nrow, ncol);
		}
		void setElev(std::vector<double> el) {
			Elev = el;
		}

		void setV(double x, double y, double z, int Nlay, int Nt) {
			vx.set(Nlay, Nt, x);
			vy.set(Nlay, Nt, y);
			vz.set(Nlay, Nt, z);
		}

		double getV(int dim, int lay, int time) {
			if (dim == 0) 
				return vx(lay, time);
			else if (dim == 1)
				return vy(lay, time);
			else if (dim == 2)
				return vz(lay, time);
			return 0;
		}

		double getElev(int lay) {
			return Elev[lay];
		}

	private:
		ic::vec3 pos;
		std::vector<double> Elev;
		mat2 vx;
		mat2 vy;
		mat2 vz;
	};

	class iwfmVel : public ICHNOS::velocityField {
	public:
		iwfmVel(boost::mpi::communicator& world_in);
		void readVelocityField(std::string vf_file, int nPnts);
		void calcVelocity(ic::vec3& vel,
                          std::vector<int>& ids,
                          std::vector<double>& weights,
                          double tm = 0);
	private:
		//ic::pointCloud<int> xyid;
		//std::unique_ptr<nano_kd_tree_int> xy_tree;
		bool readXYid(std::string prefix, std::string suffix);
		bool readElevation(std::string prefix, std::string suffix);
		bool readVelocities(std::string prefix, std::string suffix);
		std::map<int, xynode> ND;
		int Nlayers = 0;
		int Ntsteps = 0;
		ICHNOS::SingletonGenerator* RG = RG->getInstance();
	};

	iwfmVel::iwfmVel(boost::mpi::communicator& world_in)
		:
		velocityField(world_in)
	{}

	void iwfmVel::readVelocityField(std::string vf_file, int nPnts) {
		std::cout << "Velocity configuration file: " << vf_file << std::endl;
		po::options_description velocityFieldOptions("Velocity field options");
		po::variables_map vm_vfo;
		velocityFieldOptions.add_options()
			("Prefix", po::value<std::string>(), "Prefix for the input files.")
			("Suffix", po::value<std::string>(), "Suffix for the input files.")
			("Nlayers", po::value<int>()->default_value(4), "Number of layers.")
			("Ntsteps", po::value<int>()->default_value(120), "Number of time steps.")
			("OwnerThreshold", po::value<double>()->default_value(0.75), "Threshold for the processor ownership")
			;


		po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions), vm_vfo);
		OwnerThreshold = vm_vfo["OwnerThreshold"].as<double>();
		InterpolateOutsideDomain = true;

		std::string prefix = vm_vfo["Prefix"].as<std::string>();
		std::string suffix;
		if (vm_vfo.count("Suffix")) {
			suffix = vm_vfo["Suffix"].as<std::string>();
		}
		else
			suffix = ".ich";
		Nlayers = vm_vfo["Nlayers"].as<int>();
		Ntsteps = vm_vfo["Ntsteps"].as<int>();

		readXYid(prefix, suffix);

		//xy_tree = std::unique_ptr<nano_kd_tree_int>(new nano_kd_tree_int(3, xyid, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
		//xy_tree->buildIndex();
		readElevation(prefix, suffix);
		readVelocities(prefix, suffix);

	}

	bool iwfmVel::readXYid(std::string prefix, std::string suffix) {
		std::cout << "\tReading Nodes file ..." << std::endl;
		std::string filename = prefix + "XY_" + std::to_string(world.rank()) + suffix;
		std::ifstream datafile(filename.c_str());
		if (!datafile.good()) {
			std::cout << "Can't open the file" << filename << std::endl;
			return false;
		}
		else {
			std::string line;
			int id, iproc;
			ic::vec3 p;
			// Read the firt line with the parameters
			getline(datafile, line);
			std::istringstream inp(line.c_str());
			//inp >> xyid.Radious;
			//inp >> xyid.Power;
			//inp >> xyid.Threshold;

			while (getline(datafile, line)) {
				std::istringstream inp1(line.c_str());
				inp1 >> id;
				inp1 >> p.x;
				inp1 >> p.y;
				inp1 >> iproc;
				//xyid.proc.push_back(iproc);
				//xyid.pts.push_back(p);
				//xyid.vel.push_back(id);
				xynode xynd(p);
				xynd.init(Nlayers, Ntsteps);
				ND.insert(std::make_pair(id, xynd));
			}
			datafile.close();
			return true;
		}
	}

	bool iwfmVel::readElevation(std::string prefix, std::string suffix) {
		std::cout << "\tReading Elevation file ..." << std::endl;
		std::string filename = prefix + "Elev_" + std::to_string(world.rank()) + suffix;
		std::ifstream datafile(filename.c_str());
		if (!datafile.good()) {
			std::cout << "Can't open the file" << filename << std::endl;
			return false;
		}
		else {
			std::string line;
			double z;
			int id;
			std::map<int, xynode>::iterator it;
			while (getline(datafile, line)) {
				std::vector<double> zv;
				std::istringstream inp(line.c_str());
				inp >> id;
				for (int lay = 0; lay < Nlayers+1; ++lay) {
					inp >> z;
					zv.push_back(z);
				}
				it = ND.find(id);
				if (it != ND.end()) {
					it->second.setElev(zv);
				}
				else {
					std::cout << "The elevation file containts a node with id " << id << " but the XY file does not" << std::endl;
					return false;
				}
			}
			datafile.close();
			return true;
		}
	}
	bool iwfmVel::readVelocities(std::string prefix, std::string suffix) {
		for (int lay = 0; lay < Nlayers; ++lay) {
			std::cout << "\tReading velocities for layer " << lay << " ..." << std::endl;
			std::string filenameVX = prefix + "VX_LAY" + std::to_string(lay) + "_" + std::to_string(world.rank()) + suffix;
			std::string filenameVY = prefix + "VY_LAY" + std::to_string(lay) + "_" + std::to_string(world.rank()) + suffix;
			std::string filenameVZ = prefix + "VZ_LAY" + std::to_string(lay) + "_" + std::to_string(world.rank()) + suffix;
			std::ifstream datafileVX(filenameVX.c_str());
			std::ifstream datafileVY(filenameVY.c_str());
			std::ifstream datafileVZ(filenameVZ.c_str());
			bool bvx = true;
			bool bvy = true;
			bool bvz = true;
			if (!datafileVX.good()) {
				std::cout << "Can't open the file" << filenameVX << std::endl;
				bvx = false;
				return false;
			}
			if (!datafileVY.good()) {
				std::cout << "Can't open the file" << filenameVY << std::endl;
				bvy = false;
				return false;
			}
			if (!datafileVZ.good()) {
				std::cout << "Can't open the file" << filenameVZ << std::endl;
				bvz = false;
				return false;
			}

			if (bvx && bvy && bvz) {
				int idx, idy, idz;
				std::string lineX, lineY, lineZ;
				double vx, vy, vz;
				std::map<int, xynode>::iterator it;
				while (getline(datafileVX, lineX)) {
					getline(datafileVY, lineY);
					getline(datafileVZ, lineZ);
					std::istringstream inpX(lineX.c_str());
					std::istringstream inpY(lineY.c_str());
					std::istringstream inpZ(lineZ.c_str());
					inpX >> idx;
					inpY >> idy;
					inpZ >> idz;
					if (idx != idy || idx != idz || idy != idz) {
						std::cout << "The velocity files are not written in the same order" << std::endl;
						return false;
					}
					it = ND.find(idx);
					if (it != ND.end()) {
						for (int i = 0; i < Ntsteps; ++i) {
							inpX >> vx;
							inpY >> vy;
							inpZ >> vz;
							it->second.setV(vx, vy, vz, lay, i);
						}
					}
					else {
						std::cout << "Unable to find " << idx << " In the map" << std::endl;
						return false;
					}
				}
			}
			datafileVX.close();
			datafileVY.close();
			datafileVZ.close();
		}
		return true;
	}

	void iwfmVel::calcVelocity(ic::vec3& vel,
                               std::vector<int>& ids,
                               std::vector<double>& weights,
                               double tm) {
		std::map<int, double> id_dst;
		//ICHNOS::interpolateIntTree(id_dst, proc_map, xy_tree, p);

//		// Find the layer that the point is in
//		std::map<int, double>::iterator it;
//		std::map<int, xynode>::iterator itnd;
//		std::vector< std::map<int, xynode>::iterator> itvec;
//		std::vector<double> weights;
//		double zLayUp = 0.0;
//		double zLayDown = 0.0;
//		double zvel, zvelup;
//
//		double w1, w2;
//		int l1, l2;
//		bool usew2 = false;
//		for (it = id_dst.begin(); it != id_dst.end(); ++it) {
//			itnd = ND.find(it->first);
//			if (itnd != ND.end()) {
//				itvec.push_back(itnd);
//				weights.push_back(it->second);
//			}
//		}
//
//		for (int lay = 0; lay < Nlayers; ++lay) {
//			for (unsigned int i = 0; i < itvec.size(); ++i) {
//				if (lay == 0)
//					zLayUp += itvec[i]->second.getElev(lay)* weights[i];
//				zLayDown += itvec[i]->second.getElev(lay+1) * weights[i];
//			}
//
//			zvel = (zLayUp + zLayDown) / 2;
//			if (lay == 0) {
//				if (p.z > zvel) {
//					l1 = 0; w1 = 1.0;
//					l2 = 1; w2 = 0.0;
//					break;
//				}
//				else {
//					zvelup = zvel;
//				}
//			}
//			else {
//				if (p.z < zvelup && p.z > zvel) {
//					l1 = lay-1;
//					w1 = (p.z - zvel) / (zvelup - zvel);
//					l2 = lay;
//					w2 = 1.0 - w1;
//					usew2 = true;
//					break;
//				}
//				else {
//					if (lay == Nlayers - 1) {
//						if (p.z < zvel) {
//							l1 = lay; w1 = 1.0;
//							l2 = 1; w2 = 0.0;
//							break;
//						}
//						else {
//							std::cout << " The interpolation of (" << p.x << "," << p.y << "," << p.z << ") does not seem to fit into the elevation datas" << std::endl;
//						}
//					}
//					else {
//						zvelup = zvel;
//					}
//				}
//			}
//			zLayUp = zLayDown;
//			zLayDown = 0.0;
//		}
//
//		//std::cout << "we landed here safely" << std::endl;
//		vel.zero();
//		int r = -9;
//		while (r < 0 || r > Ntsteps - 1)
//			r = RG->randomNumber(0, Ntsteps);
//		for (unsigned int i = 0; i < itvec.size(); ++i) {
//			vel.x += itvec[i]->second.getV(0, l1, r) * weights[i] * w1;
//			vel.y += itvec[i]->second.getV(1, l1, r) * weights[i] * w1;
//			vel.z += itvec[i]->second.getV(2, l1, r) * weights[i] * w1;
//			if (usew2) {
//				vel.x += itvec[i]->second.getV(0, l2, r) * weights[i] * w2;
//				vel.y += itvec[i]->second.getV(1, l2, r) * weights[i] * w2;
//				vel.z += itvec[i]->second.getV(2, l2, r) * weights[i] * w2;
//			}
//		}
	}
}