#pragma once

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include <nanoflann.hpp>

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
		void readVelocityField(std::string vf_file);
		void calcVelocity(ic::vec3& vel, std::map<int, double>& proc_map, ic::vec3& p);
	private:
		ic::pointCloud<int> xyid;
		std::unique_ptr<nano_kd_tree_int> xy_tree;
		void readXYid(std::string prefix, std::string suffix);
		void readElevation(std::string prefix, std::string suffix);
		bool readVelocities(std::string prefix, std::string suffix);
		std::map<int, xynode> ND;
		int Nlayers;
		int Ntsteps;
	};

	iwfmVel::iwfmVel(boost::mpi::communicator& world_in)
		:
		velocityField(world_in)
	{}

	void iwfmVel::readVelocityField(std::string vf_file) {
		po::options_description velocityFieldOptions("Velocity field options");
		po::variables_map vm_vfo;
		velocityFieldOptions.add_options()
			("Prefix", po::value<std::string>(), "Prefix for the input files.")
			("Suffix", po::value<std::string>(), "Suffix for the input files.")
			("Nlayers", po::value<std::string>(), "Number of layers.")
			("Ntsteps", po::value<std::string>(), "Number of time steps.")
			("OwnerThreshold", po::value<double>()->default_value(0.75), "Threshold for the processor ownership")
			;

		po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions), vm_vfo);
		OwnerThreshold = vm_vfo["OwnerThreshold"].as<double>();
		InterpolateOutsideDomain = true;

		std::string prefix = vm_vfo["pcPrefix"].as<std::string>();
		std::string suffix;
		if (vm_vfo.count("Prefix")) {
			suffix = vm_vfo["Suffix"].as<std::string>();
		}
		else
			suffix = ".ich";

		readXYid(prefix, suffix);

		xy_tree = std::unique_ptr<nano_kd_tree_int>(new nano_kd_tree_int(3, xyid, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
		xy_tree->buildIndex();
		readElevation(prefix, suffix);
		readVelocities(prefix, suffix);

	}

	void iwfmVel::readXYid(std::string prefix, std::string suffix) {
		std::string filename = prefix + "XY_" + std::to_string(world.rank()) + suffix;
		std::ifstream datafile(filename.c_str());
		if (!datafile.good()) {
			std::cout << "Can't open the file" << filename << std::endl;
		}
		else {
			std::string line;
			int id, iproc;
			ic::vec3 p;
			while (getline(datafile, line)) {
				std::istringstream inp(line.c_str());
				inp >> id;
				inp >> p.x;
				inp >> p.y;
				inp >> iproc;
				xyid.proc.push_back(iproc);
				xyid.pts.push_back(p);
				xyid.vel.push_back(id);
				xynode xynd(p);
				xynd.init(Nlayers, Ntsteps);
				ND.insert(std::make_pair(id, xynd));
			}
			datafile.close();
		}
	}

	void iwfmVel::readElevation(std::string prefix, std::string suffix) {
		std::string filename = prefix + "Elev_" + std::to_string(world.rank()) + suffix;
		std::ifstream datafile(filename.c_str());
		if (!datafile.good()) {
			std::cout << "Can't open the file" << filename << std::endl;
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
				for (int lay = 0; lay < Nlayers; ++lay) {
					inp >> z;
					zv.push_back(z);
				}
				it = ND.find(id);
				if (it != ND.end()) {
					it->second.setElev(zv);
				}
			}
			datafile.close();
		}
	}
	bool iwfmVel::readVelocities(std::string prefix, std::string suffix) {
		for (int lay = 0; lay < Nlayers; ++lay) {
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
						std::cout << "THe velocity files are not written in the same order" << std::endl;
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
	}

	void iwfmVel::calcVelocity(ic::vec3& vel, std::map<int, double>& proc_map, ic::vec3& p) {

	}
}