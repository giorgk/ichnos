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
		int proc = -9;
		//double dist;
		ic::vec3 p;
		ic::vec3 vel;
	};

	class npsatVel : public ic::velocityField {
	public:
		npsatVel(boost::mpi::communicator& world_in);
		void readVelocityField(std::string vf_file);
		void calcVelocity(ic::vec3& vel, std::map<int, double>& proc_map, ic::vec3& p, double& step);
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
		double max_search_time = 0.0;
		double calc_time = 0.0;
		double max_calc_time = 0.0;
		double calc_time1 = 0.0;
		double calc_time2 = 0.0;
		int cout_times = 0;
		int FrequencyStat;
		int nids_size = 0;
		int method;
		//double max_distance = 0;
		void PrintStat();

		void add_interpolationVertex(std::map<int, npsatVeldata>& nids, ic::vec3& lp, ic::vec3& up, ic::vertex_handle& vh, ic::vec3& p, npsatVeldata& tmp);
		void calculate_BB(ic::vertex_handle& vh, ic::vec3& l, ic::vec3& u);
		double calculate_step(ic::vec3& p, ic::vec3& v, ic::vec3& l, ic::vec3& u);
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
			("Method", po::value<int>()->default_value(1), "Method to identify nearby vertices")
			;

		po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions), vm_vfo);

		OwnerThreshold = vm_vfo["OwnerThreshold"].as<double>();
		int leadZeros = vm_vfo["LeadingZeros"].as<int>();
		multiplier = vm_vfo["VelocityMultiplier"].as<double>();
		Power = vm_vfo["Power"].as<double>();
		Scale = vm_vfo["Scale"].as<double>();
		Threshold = vm_vfo["Threshold"].as<double>();
		FrequencyStat = vm_vfo["FrequencyStat"].as<int>();
		method = vm_vfo["Method"].as<int>();

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



	void npsatVel::calcVelocity(ic::vec3& vel, std::map<int, double>& proc_map, ic::vec3& p, double& step) {
		npsatVeldata tmp;
		ic::vertex_handle vh;
		ic::cell_handle ch;
		//ic::NPSAT_data nd;
		std::map<int, npsatVeldata> nids; // nearest point ids
		std::map<int, npsatVeldata>::iterator it;
		std::vector<ic::vertex_handle> inc_v; // container for the incident vertices
		ic::vec3 lp(999999999, 999999999, 999999999);
		ic::vec3 up(-999999999, -999999999, -999999999);
		ic::vec3 lbb(999999999, 999999999, 999999999);
		ic::vec3 ubb(-999999999, -999999999, -999999999);
		//max_distance = 0.0;

		auto start = std::chrono::high_resolution_clock::now();
		if (method == 1){
			vh = T.nearest_vertex(ic::cgal_point(p.x, p.y, p.z), cellhandle);
			auto finish = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = finish - start;
			search_time += elapsed.count();
			if (elapsed.count() > max_search_time)
				max_search_time = elapsed.count();

			start = std::chrono::high_resolution_clock::now();
			add_interpolationVertex(nids, lp, up, vh, p, tmp);

			//------- Test using only the incident vertices to this vertex
			if (T.is_infinite(vh))
				std::cout << "The nearest vertex is infinite" << std::endl;

			T.finite_adjacent_vertices(vh, std::back_inserter(inc_v));
			for (unsigned int ii = 0; ii < inc_v.size(); ++ii) {
				add_interpolationVertex(nids, lp, up, inc_v[ii], p, tmp);
			}
			ch = vh->cell();
			if (!T.is_infinite(ch)) {
				cellhandle = ch;
			}
		}
		else if (method == 2){
			vh = T.nearest_vertex(ic::cgal_point(p.x, p.y, p.z), cellhandle);
			auto finish = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = finish - start;
			search_time += elapsed.count();
			if (elapsed.count() > max_search_time)
				max_search_time = elapsed.count();

			start = std::chrono::high_resolution_clock::now();
			add_interpolationVertex(nids, lp, up, vh, p, tmp);

			if (T.is_infinite(vh))
				std::cout << "The nearest vertex is infinite" << std::endl;

			ch = vh->cell();
			if (!T.is_infinite(ch)) {
				cellhandle = ch;
			}

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
		}
		else if (method == 3){
			ch = T.locate(ic::cgal_point(p.x, p.y, p.z), cellhandle);
			auto finish = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = finish - start;
			search_time += elapsed.count();
			if (elapsed.count() > max_search_time)
				max_search_time = elapsed.count();
			
			if (!T.is_infinite(ch)) {
				cellhandle = ch;
			}
			for (int i = 0; i < 4; ++i){
				vh = ch->vertex(i);
				if (T.is_infinite(vh))
					continue;
				if (vh->info().id < 0)
					continue;
				add_interpolationVertex(nids, lp, up, vh, p, tmp);
			}
		}
		else if (method == 4){
			ch = T.locate(ic::cgal_point(p.x, p.y, p.z), cellhandle);
			auto finish = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = finish - start;
			search_time += elapsed.count();
			if (elapsed.count() > max_search_time)
				max_search_time = elapsed.count();
			
			if (!T.is_infinite(ch)) {
				cellhandle = ch;
			}

			for (int i = 0; i < 4; ++i){
				vh = ch->vertex(i);
				if (T.is_infinite(vh))
					continue;
				if (vh->info().id < 0 || vh->info().proc < 0)
					continue;

				calculate_BB(vh, lbb, ubb);

				add_interpolationVertex(nids, lp, up, vh, p, tmp);
				inc_v.clear();
				T.finite_adjacent_vertices(vh, std::back_inserter(inc_v));
				for (unsigned int ii = 0; ii < inc_v.size(); ++ii) {
					add_interpolationVertex(nids, lp, up, inc_v[ii], p, tmp);
				}
			}
		}
		
		if (nids.size() == 0){
			std::cout << "The number of Nodes is ZERO" << std::endl;
			vel = ic::vec3();
		}
		//auto finish1 = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> elapsed1 = finish1 - start;
		//calc_time1 += elapsed1.count();

		//auto start1 = std::chrono::high_resolution_clock::now();
		double horlen = std::min(up.x - lp.x, up.y - lp.y);
		double verlen = up.z - lp.z;

		double ratio;
		if (verlen < 0.1){
			ratio = 1.0;
		}
		else{
			ratio = horlen / verlen;
		} 

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
					if (std::isnan(w)){
						std::cout << "w is nan" << std::endl;
					}
					itd->second += w;
				}
				sumW += w;
				sumWVal = sumWVal + it->second.vel * w;
			}
		}
		if (!std::isfinite(sumWVal.x)){
			std::cout << "Nan Velocity" << std::endl;
		}

		//std::cout << std::endl;
		//finish1 = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> elapsed2 = finish1 - start1;
		//calc_time2 += elapsed2.count();
		itd = proc_map.begin();
		for (; itd != proc_map.end(); ++itd) {
			itd->second = itd->second / sumW;
		}
		if (calc_average)
			vel = sumWVal * (1 / sumW);
		if (!std::isfinite(std::abs(vel.x)))
			std::cout << vel.x << "," << vel.y << "," << vel.z << std::endl;
		//ic::interpolateVectorTree(vel, proc_map, VelocityTree, p);
		double porosity = 1.0;
		if (porType == ic::interpType::CLOUD) {
			//porosity = ic::interpolateScalarTree(PorosityTree, p);
		}
		else if (porType == ic::interpType::SCALAR)
			porosity = porosityValue;
		vel = vel * (1/porosity);

		if (method == 4){
			step = calculate_step(p, vel, lbb, ubb);
		}

		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		calc_time += elapsed.count();
		if (elapsed.count() > max_calc_time)
			max_calc_time = elapsed.count();
		cout_times++;
		nids_size = nids.size();
		PrintStat();
	}

	void npsatVel::add_interpolationVertex(std::map<int, npsatVeldata>& nids, ic::vec3& lp, ic::vec3& up, ic::vertex_handle& vh, ic::vec3& p, npsatVeldata& tmp) {
		ic::NPSAT_data nd = vh->info();
		if (nd.id < 0 || nd.proc < 0)
			return;

		if (T.is_infinite(vh))
			return;

		if (nids.find(nd.id) != nids.end())
			return;

		tmp.proc = nd.proc;
		tmp.vel = nd.v;
		ic::cgal_point tmp_p = vh->point();
		//std::cout << tmp_p << std::endl;
		tmp.p = ic::vec3(p.x - tmp_p.x(), p.y - tmp_p.y(), p.z - tmp_p.z());
		//double ddd = tmp.p.len();
		//if (ddd > max_distance){
		//	max_distance = ddd;
		//}
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

	void npsatVel::calculate_BB(ic::vertex_handle& vh, ic::vec3& l, ic::vec3& u){
		ic::cgal_point tmp_p = vh->point();
		if (tmp_p.x() > u.x)
			u.x = tmp_p.x();
		if (tmp_p.y() > u.y)
			u.y = tmp_p.y();
		if (tmp_p.z() > u.z)
			u.z = tmp_p.z();
		if (tmp_p.x() < l.x)
			l.x = tmp_p.x();
		if (tmp_p.y() < l.y)
			l.y = tmp_p.y();
		if (tmp_p.z() < l.z)
			l.z = tmp_p.z();
	}

	double npsatVel::calculate_step(ic::vec3& p, ic::vec3& v, ic::vec3& l, ic::vec3& u){
		ic::vec3 p1 = p + v * 0.1;
		double t_max = 1000000000;
		double t_min = -1000000000;
		double tmp_min, tmp_max;

		if (std::abs(v.x) < 0.000001){
			tmp_min = -1000000000;
			tmp_max =  1000000000;
		}
		else if (v.x > 0.0 ){
			tmp_min = (l.x - p.x)/(p1.x - p.x);
			tmp_max = (u.x - p.x)/(p1.x - p.x);
		}
		else if (v.x < 0.0){
			tmp_min = (u.x - p.x)/(p1.x - p.x);
			tmp_max = (l.x - p.x)/(p1.x - p.x);
		}
		if (tmp_max < t_max)
			t_max = tmp_max;
		if (tmp_min > t_min)
			t_min = tmp_min;

		if (std::abs(v.y) < 0.000001){
			tmp_min = -1000000000;
			tmp_max =  1000000000;
		}
		else if (v.y > 0.0 ){
			tmp_min = (l.y - p.y)/(p1.y - p.y);
			tmp_max = (u.y - p.y)/(p1.y - p.y);
		}
		else if (v.y < 0.0){
			tmp_min = (u.y - p.y)/(p1.y - p.y);
			tmp_max = (l.y - p.y)/(p1.y - p.y);
		}
		if (tmp_max < t_max)
			t_max = tmp_max;
		if (tmp_min > t_min)
			t_min = tmp_min;

		if (std::abs(v.z) < 0.000001){
			tmp_min = -1000000000;
			tmp_max =  1000000000;
		}
		else if (v.z > 0.0 ){
			tmp_min = (l.z - p.z)/(p1.z - p.z);
			tmp_max = (u.z - p.z)/(p1.z - p.z);
		}
		else if (v.z < 0.0){
			tmp_min = (u.z - p.z)/(p1.z - p.z);
			tmp_max = (l.z - p.z)/(p1.z - p.z);
		}
		if (tmp_max < t_max)
			t_max = tmp_max;
		if (tmp_min > t_min)
			t_min = tmp_min;

		ic::vec3 pmin = p*(1-t_min) + p1 * t_min;
		ic::vec3 pmax = p*(1-t_max) + p1 * t_max;
		return pmin.distance(pmax.x, pmax.y, pmax.z)/2;
	}

	void npsatVel::PrintStat() {
		if (cout_times > FrequencyStat) {
			std::cout << "Search time: " << std::fixed << std::setprecision(15) << search_time/static_cast<double>(cout_times) << ", ("; // std::endl;
			std::cout << max_search_time << "), ";
			std::cout << "Velocity Calc time: " << std::fixed << std::setprecision(15) << calc_time/static_cast<double>(cout_times) <<  ", (";// std::endl;
			std::cout << max_calc_time << "), ";
			//std::cout << std::endl << std::flush;
			//std::cout << "max: " << max_distance << ", ";
			std::cout << "Number of nodes: " << nids_size << std::endl;
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

