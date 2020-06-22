#pragma once

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
//#include <nanoflann.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include "ichnos_structures.h"
#include "velocity_base.h"

#include <chrono>
#include <ctime>

namespace po = boost::program_options;
namespace ic = ICHNOS;


namespace NPSAT {

	struct npsatVeldata {
		int proc;
		//double dist;
		ic::vec3 p;
		ic::vec3 vel;
	};

	class npsatVel : public ic::velocityField {
	public:
		npsatVel(boost::mpi::communicator& world_in);
		void readVelocityField(std::string vf_file);
		void calcVelocity(ic::vec3& vel, std::map<int, double>& proc_map, ic::vec3& p);
	private:
		ic::cgal_Delaunay T;
		ic::cell_handle cellhandle;
		bool bCellKnown = false;
		//ic::pointCloud<ic::vec3> VelocityCloud;
		//std::unique_ptr<nano_kd_tree_vector> VelocityTree;
		ic::interpType porType;
		//ic::pointCloud<double> PorosityCloud;
		double porosityValue = 1.0;
		//std::unique_ptr<nano_kd_tree_scalar> PorosityTree;
		double multiplier = 1.0;
		double Power;
		double Scale = 1.0;
		double Threshold;
		double search_time = 0.0;
		double calc_time = 0.0;
		double calc_time1 = 0.0;
		double calc_time2 = 0.0;
		int cout_times = 0;
		int FrequencyStat;
		int nids_size = 0;
		void PrintStat();

		void add_interpolationVertex(std::map<int, npsatVeldata>& nids, ic::vec3& lp, ic::vec3& up, ic::vertex_handle& vh, ic::vec3& p, npsatVeldata& tmp);
	};

	npsatVel::npsatVel(boost::mpi::communicator& world_in)
		:
		velocityField(world_in)
	{
		InterpolateOutsideDomain = true;
	}

	void npsatVel::readVelocityField(std::string vf_file) {
		if (world.rank() == 0)
			std::cout << "Reading data..." << std::endl;
		po::options_description velocityFieldOptions("Velocity field options");
		po::variables_map vm_vfo;
		velocityFieldOptions.add_options()
			("Prefix", po::value<std::string>(), "Prefix for the filename") 
			("LeadingZeros", po::value<int>()->default_value(4), "e.g 0002->4, 000->3")
			("Suffix", po::value<std::string>(), "ending of file after procid")
			("OwnerThreshold", po::value<double>()->default_value(0.75), "Threshold for the processor ownership")
			("VelocityMultiplier", po::value<double>()->default_value(1), "This is a multiplier to scale velocity")
			("Porosity", po::value<std::string>(), "Porocity. Either a file or a single number")
			("Power", po::value<double>()->default_value(3.0), "Power of the IDW interpolation")
			("Scale", po::value<double>()->default_value(1.0), "Scale the domain before velocity calculation")
			("Threshold", po::value<double>()->default_value(0.001), "Threshold of distance of IDW")
			("FrequencyStat", po::value<int>()->default_value(20), "Frequency of printing stats")
			;

		po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions), vm_vfo);

		OwnerThreshold = vm_vfo["OwnerThreshold"].as<double>();
		int leadZeros = vm_vfo["LeadingZeros"].as<int>();
		multiplier = vm_vfo["VelocityMultiplier"].as<double>();
		Power = vm_vfo["Power"].as<double>();
		Scale = vm_vfo["Scale"].as<double>();
		Threshold = vm_vfo["Threshold"].as<double>();
		FrequencyStat = vm_vfo["FrequencyStat"].as<int>();

		std::string prefix = vm_vfo["Prefix"].as<std::string>();
		std::string suffix;
		if (vm_vfo.count("Suffix")) {
			suffix = vm_vfo["Suffix"].as<std::string>();
		}
		else
			suffix = ".ich";

		std::string filename = prefix + ic::num2Padstr(world.rank(), leadZeros) + suffix;
		std::vector<std::pair<ic::cgal_point, ic::NPSAT_data>> npsat_data;
		ic::READ::readNPSATVelocity(filename, npsat_data, multiplier);
		auto start = std::chrono::high_resolution_clock::now();
		T.insert(npsat_data.begin(), npsat_data.end());
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Triangulation Building time: " << elapsed.count() << std::endl;
		std::cout << "Number of vertices: " << T.number_of_vertices() << std::endl;
		std::cout << "Number of cells: " << T.number_of_cells() << std::endl;


		//ic::READ::readVelocityFieldFile(filename, VelocityCloud);
		//this->InterpolateOutsideDomain = VelocityCloud.InterpolateOutside;
		//this->VelocityTree = std::unique_ptr<nano_kd_tree_vector>(new nano_kd_tree_vector(3, VelocityCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
		//VelocityTree->buildIndex();

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
					//porType = ic::interpType::CLOUD;
					//ic::READ::read2DscalarField(porfile, PorosityCloud);
					//this->PorosityTree = std::unique_ptr<nano_kd_tree_scalar>(new nano_kd_tree_scalar(3, PorosityCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
					//PorosityTree->buildIndex();
				}
			}
		}
		else {
			porType = ic::interpType::INGORE;
		}
	}



	void npsatVel::calcVelocity(ic::vec3& vel, std::map<int, double>& proc_map, ic::vec3& p) {
		npsatVeldata tmp;
		ic::vertex_handle vh;
		ic::cell_handle ch;
		//ic::NPSAT_data nd;
		std::map<int, npsatVeldata> nids; // nearest point ids
		std::map<int, npsatVeldata>::iterator it;
		std::vector<ic::vertex_handle> inc_v; // container for the incident vertices
		ic::vec3 lp(999999999, 999999999, 999999999);
		ic::vec3 up(-999999999, -999999999, -999999999);

		auto start = std::chrono::high_resolution_clock::now();
		vh = T.nearest_vertex(ic::cgal_point(p.x, p.y, p.z), cellhandle);
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		search_time += elapsed.count();

		start = std::chrono::high_resolution_clock::now();
		add_interpolationVertex(nids, lp, up, vh, p, tmp);

		ch = vh->cell();
		//std::cout << T.is_infinite(ch) << " " << std::flush;
		//cellhandle = T.locate(ic::cgal_point(p.x, p.y, p.z), cellhandle);
		//bool debug_this = false;
		if (!T.is_infinite(ch)) {
			cellhandle = ch;
		}
		//else{
		//	debug_this = true;
		//}
		for (int i = 0; i < 4; ++i) {
			vh = ch->vertex(i);
			if (T.is_infinite(vh))
				continue;
			//if (debug_this){
			//	std::cout << vh->info().id << std::endl;
			//	std::cout << T.is_infinite(vh) << std::endl;
			//}
			if (vh->info().id < 0)
				continue;
			
			add_interpolationVertex(nids, lp, up, vh, p, tmp);
			inc_v.clear();
			T.finite_adjacent_vertices(vh, std::back_inserter(inc_v));
			for (unsigned int ii = 0; ii < inc_v.size(); ++ii) {
				add_interpolationVertex(nids, lp, up, inc_v[ii], p, tmp);
			}
		}
		auto finish1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed1 = finish1 - start;
		calc_time1 += elapsed1.count();

		auto start1 = std::chrono::high_resolution_clock::now();
		double horlen = std::min(up.x - lp.x, up.y - lp.y);
		double verlen = up.z - lp.z;

		double ratio = horlen / verlen;

		double dist, w;
		double sumW = 0;
		ic::vec3 sumWVal;
		std::map<int, double>::iterator itd;
		bool calc_average = true;
		for (it = nids.begin(); it != nids.end(); ++it) {
			//std::cout << it->first << " ";
			it->second.p.z = it->second.p.z * ratio * Scale + it->second.p.z * (1.0 - Scale);
			dist = std::sqrt(it->second.p.x * it->second.p.x + it->second.p.y * it->second.p.y + it->second.p.z * it->second.p.z);

			if (dist < Threshold) {
				vel = it->second.vel;
				calc_average = false;
			}
			else {
				w = 1 / std::pow(dist, Power);
				itd = proc_map.find(it->second.proc);
				if (itd == proc_map.end()) {
					proc_map.insert(std::pair<int, double>(it->second.proc, w));
				}
				else {
					itd->second += w;
				}
				sumW += w;
				sumWVal = sumWVal + it->second.vel * w;
			}
		}
		//std::cout << std::endl;
		finish1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed2 = finish1 - start1;
		calc_time2 += elapsed2.count();
		itd = proc_map.begin();
		for (; itd != proc_map.end(); ++itd) {
			itd->second = itd->second / sumW;
		}
		if (calc_average)
			vel = sumWVal * (1 / sumW);

		//ic::interpolateVectorTree(vel, proc_map, VelocityTree, p);
		double porosity = 1.0;
		if (porType == ic::interpType::CLOUD) {
			//porosity = ic::interpolateScalarTree(PorosityTree, p);
		}
		else if (porType == ic::interpType::SCALAR)
			porosity = porosityValue;
		vel = vel * (1/porosity);

		finish = std::chrono::high_resolution_clock::now();
		elapsed = finish - start;
		calc_time += elapsed.count();
		cout_times++;
		nids_size = nids.size();
		PrintStat();
	}

	void npsatVel::add_interpolationVertex(std::map<int, npsatVeldata>& nids, ic::vec3& lp, ic::vec3& up, ic::vertex_handle& vh, ic::vec3& p, npsatVeldata& tmp) {
		ic::NPSAT_data nd = vh->info();
		if (nd.id < 0)
			return;

		if (nids.find(nd.id) != nids.end())
			return;

		tmp.proc = nd.proc;
		tmp.vel = nd.v;
		ic::cgal_point tmp_p = vh->point();
		//std::cout << tmp_p << std::endl;
		tmp.p = ic::vec3(p.x - tmp_p.x(), p.y - tmp_p.y(), p.z - tmp_p.z());
		//tmp.dist = p.distance(tmp_p.x(), tmp_p.y(), tmp_p.z());
		nids.insert(std::make_pair(nd.id, tmp));
		//nids[nd.id] = tmp;
		if (tmp_p.x() > up.x)
			up.x = tmp_p.x();
		if (tmp_p.y() > up.y)
			up.y = tmp_p.y();
		if (tmp_p.z() > up.z)
			up.z = tmp_p.z();
		if (tmp_p.x() < lp.x)
			lp.x = tmp_p.x();
		if (tmp_p.y() < lp.y)
			lp.y = tmp_p.y();
		if (tmp_p.z() < lp.z)
			lp.z = tmp_p.z();
	}

	void npsatVel::PrintStat() {
		if (cout_times > FrequencyStat) {
			std::cout << "Search time: " << std::fixed << std::setprecision(15) << search_time/static_cast<double>(cout_times) << ", "; // std::endl;
			std::cout << "Velocity Calc time: " << std::fixed << std::setprecision(15) << calc_time/static_cast<double>(cout_times) <<  ", ";// std::endl;
			std::cout << std::endl << std::flush;
			//std::cout << "Number of nodes: " << nids_size << std::endl;
			//std::cout << "Velocity Calc time1: " << std::fixed << std::setprecision(15) << calc_time1 / static_cast<double>(cout_times) << std::endl;
			//std::cout << "Velocity Calc time2: " << std::fixed << std::setprecision(15) << calc_time2 / static_cast<double>(cout_times) << std::endl;
 			cout_times = 0;
			search_time = 0.0;
			calc_time = 0.0;
			calc_time1 = 0.0;
			calc_time2 = 0.0;
		}
	}
}

