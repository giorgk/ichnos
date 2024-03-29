#pragma once

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include "ichnos_structures.h"
#include "velocity_base.h"

#include <chrono>
#include <ctime>
#include <utility>

namespace po = boost::program_options;
namespace ic = ICHNOS;


namespace NPSAT {

	/*struct npsatVeldata {
		int proc = -9;
		//double dist;
		ic::vec3 p;
		ic::vec3 vel;
	};*/

	class npsatVel : public ic::velocityField {
	public:
		npsatVel(boost::mpi::communicator& world_in);
		bool readVelocityField(std::string vf_file);
		void calcVelocity(ic::vec3& vel,
                          std::vector<int>& ids,
                          std::vector<double>& weights,
                          double tm = 0);
		void reset();
		void updateStep(double& step);
	private:

		//ic::search_tree Tree;

		//ic::cgal_Delaunay T;
		//ic::cell_handle cellhandle;
		//bool bCellKnown = false;
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
		//double search_time = 0.0;
		//double max_search_time = 0.0;
		double calc_time = 0.0;
		double max_calc_time = 0.0;
		//double calc_time1 = 0.0;
		//double calc_time2 = 0.0;
		int count_times = 0;
		int FrequencyStat;
		int nids_size = 0;
		//int method;
		//double max_distance = 0;

		bool bIsInitialized = false;
		double initial_diameter = 640;
		double initial_ratio = 20;
		double diameter;
		double ratio;
		double search_mult = 2.5;
		double n_steps = 4.0;
		ic::vec3 ll,uu,pp,vv;

		//void add_interpolationVertex(std::map<int, npsatVeldata>& nids, ic::vec3& lp, ic::vec3& up, ic::vertex_handle& vh, ic::vec3& p, npsatVeldata& tmp);
		void calculate_BB(ic::vertex_handle& vh, ic::vec3& l, ic::vec3& u);
		double calculate_step(ic::vec3& p, ic::vec3& v, ic::vec3& l, ic::vec3& u);
	};

	npsatVel::npsatVel(boost::mpi::communicator& world_in)
		:
		velocityField(world_in)
	{
		InterpolateOutsideDomain = true;
	}

	void npsatVel::reset(){
		bIsInitialized = false;
		diameter = initial_diameter;
		ratio = initial_ratio;
	}

	bool npsatVel::readVelocityField(std::string vf_file) {
		if (world.rank() == 0)
			std::cout << "Reading data..." << std::endl;
		po::options_description velocityFieldOptions("Velocity field options");
		po::variables_map vm_vfo;
		velocityFieldOptions.add_options()
			// Velocity parameters
			("Velocity.Prefix", po::value<std::string>(), "Prefix for the filename") 
			("Velocity.LeadingZeros", po::value<int>()->default_value(4), "e.g 0002->4, 000->3")
			("Velocity.Suffix", po::value<std::string>(), "ending of file after procid")
			("Velocity.Multiplier", po::value<double>()->default_value(1), "This is a multiplier to scale velocity")
			("Velocity.Scale", po::value<double>()->default_value(1.0), "Scale the domain before velocity calculation")
			("Velocity.Power", po::value<double>()->default_value(3.0), "Power of the IDW interpolation")
			("Velocity.InitDiameter", po::value<double>()->default_value(5000), "Initial diameter")
			("Velocity.InitRatio", po::value<double>()->default_value(1), "Initial ratio")

			// Porosity parameters
			("Porosity.Value", po::value<std::string>(), "Porosity. Either a file or a single number")

			//General
			("General.OwnerThreshold", po::value<double>()->default_value(0.75), "Threshold for the processor ownership")
			("General.Threshold", po::value<double>()->default_value(0.001), "Threshold of distance of IDW")
			("General.FrequencyStat", po::value<int>()->default_value(20), "Frequency of printing stats")
			("General.Nsteps", po::value<double>()->default_value(4), "Number of steps")
			;

		po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions), vm_vfo);

		std::vector<ic::cgal_point_3> pp;
		std::vector<ic::NPSAT_data> dd;

		{// General
			OwnerThreshold = vm_vfo["General.OwnerThreshold"].as<double>();
			Threshold = vm_vfo["General.Threshold"].as<double>();
			FrequencyStat = vm_vfo["General.FrequencyStat"].as<int>();
			n_steps = vm_vfo["General.Nsteps"].as<double>();
		}


		{// Velocity parameters
			multiplier = vm_vfo["Velocity.Multiplier"].as<double>();
			Scale = vm_vfo["Velocity.Scale"].as<double>();
			Power = vm_vfo["Velocity.Power"].as<double>();
			initial_diameter = vm_vfo["Velocity.InitDiameter"].as<double>();
			initial_ratio = vm_vfo["Velocity.InitRatio"].as<double>();

			std::string prefix = vm_vfo["Velocity.Prefix"].as<std::string>();
			std::string suffix;
			if (vm_vfo.count("Velocity.Suffix")) {
				suffix = vm_vfo["Velocity.Suffix"].as<std::string>();
			}
			else {
                suffix = ".ich";
            }
			int leadZeros = vm_vfo["Velocity.LeadingZeros"].as<int>();
			std::string filename = prefix + ic::num2Padstr(/*dbg_rank*/world.rank(), leadZeros) + suffix;
			ic::READ::readNPSATVelocity(filename, pp, dd, multiplier);
		}

		{//Build tree
			auto start = std::chrono::high_resolution_clock::now();
			//Tree.insert(boost::make_zip_iterator(boost::make_tuple( pp.begin(),dd.begin() )),
    		//			boost::make_zip_iterator(boost::make_tuple( pp.end(),dd.end() ) )  );
			//Tree.build();
			auto finish = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = finish - start;
			std::cout << "Point Set Building time: " << elapsed.count() << std::endl;

		}

		if (vm_vfo.count("Porosity.Value")) {
			std::string porfile = vm_vfo["Porosity.Value"].as<std::string>();
			if (porfile.empty()) {
				porType = ic::interpType::INGORE;
			}
			else {
				if (ic::is_input_scalar(porfile)) {
					porType = ic::interpType::SCALAR;
					porosityValue = std::stod(porfile);
				}
				else {
				    // TODO
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
		return true;
	}



	void npsatVel::calcVelocity(ic::vec3& vel,
                                std::vector<int>& ids,
                                std::vector<double>& weights,
                                double time) {
	    ic::vec3 p;
		// If this is the first point of this streamline we will carry out one additional range search
		ll.zero();
		uu.zero();
		pp.zero();
		vv.zero();
		if (!bIsInitialized){
			ic::calculate_search_box(p,ll,uu,diameter,ratio,search_mult);
			ic::cgal_point_3 llp(ll.x, ll.y, ll.z);
			ic::cgal_point_3 uup(uu.x, uu.y, uu.z);
			ic::Fuzzy_iso_box fib(llp,uup, 0.0);
			std::vector<boost::tuples::tuple<ic::cgal_point_3, ic::NPSAT_data>> tmp;
			//std::back_insert_iterator<std::vector<std::tuple<ic::cgal_point_3, ic::NPSAT_data>>> tmp_it(tmp);
			//Tree.search(std::back_inserter(tmp), fib);
			if (tmp.size() == 0){
				vel = ic::vec3(-99999,-99999,-99999);
				return;
			}
			// Find the closest point
			double mindist = 1000000000;
			double tmp_diam;
			double tmp_ratio;
			ic::NPSAT_data closest_point_data;
			for (unsigned int i = 0; i < tmp.size(); ++i){
				double dist = p.distance(tmp[i].get<0>().x(), tmp[i].get<0>().y(), tmp[i].get<0>().z());
				if (dist < mindist){
					mindist = dist;
					tmp_diam = tmp[i].get<1>().diameter;
					tmp_ratio = tmp[i].get<1>().ratio;
					//std::cout << tmp[1].get<1>().id << std::endl;
				}
			}
			diameter = tmp_diam;
			ratio = tmp_ratio;
			bIsInitialized = true;
		}

		{
			auto start = std::chrono::high_resolution_clock::now();
			std::map<int, double>::iterator itd;
			std::vector<boost::tuples::tuple<ic::cgal_point_3, ic::NPSAT_data>> tmp;
			while (true){
				tmp.clear();
				ic::calculate_search_box(p,ll,uu,diameter,ratio,search_mult);
				ic::cgal_point_3 llp(ll.x, ll.y, ll.z);
				ic::cgal_point_3 uup(uu.x, uu.y, uu.z);
				ic::Fuzzy_iso_box fib(llp,uup, 0.0);
				//Tree.search(std::back_inserter(tmp), fib);
				if (tmp.size() >= 3){
					break;
				}
				else{
					diameter = diameter*1.5;
					if (diameter > initial_diameter){
						vel = ic::vec3(-99999,-99999,-99999);
						return;
					}
				}
			}

			// por_xyz refers to point translated to origin
			double porx, pory, porz, scaled_dist, actual_dist, w;
			bool calc_average = true;
			double sumW = 0;
			ic::vec3 sumWVal;
			double mindist = 1000000000;
			double tmp_diam;
			double tmp_ratio;

			for (unsigned int i = 0; i < tmp.size(); ++i){
				porx = p.x - tmp[i].get<0>().x();
				pory = p.y - tmp[i].get<0>().y();
				porz = p.z - tmp[i].get<0>().z();
				actual_dist = std::sqrt(porx * porx + pory * pory + porz * porz);
				if (actual_dist < mindist){
					mindist = actual_dist;
					tmp_diam = tmp[i].get<1>().diameter;
					tmp_ratio = tmp[i].get<1>().ratio;
				}

				porz = porz * ratio * Scale + porz*(1 - Scale);

				scaled_dist = std::sqrt(porx * porx + pory * pory + porz * porz);


				if (actual_dist < Threshold){
					vel = tmp[i].get<1>().v;
					calc_average = false;
				}
				else{
					w = 1 / std::pow(scaled_dist, Power);
//					itd = proc_map.find(tmp[i].get<1>().proc);
//					if (itd == proc_map.end()) {
//						proc_map.insert(std::pair<int, double>(tmp[i].get<1>().proc, w));
//					}
//					else{
//						itd->second += w;
//					}
					sumW += w;
					sumWVal = sumWVal + tmp[i].get<1>().v*w;
				}
			}
			if (tmp.size() > 50){
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
			calc_time += elapsed.count();
			if (elapsed.count() > max_calc_time)
			max_calc_time = elapsed.count();
			//std::cout << "Range search time: " << std::fixed << std::setprecision(15) << elapsed.count() << " N points: " << tmp.size() <<  std::endl;
		}


/*
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
			std::cout << "Search time Tria: " << elapsed.count() << std::endl;
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

		double tmp_ratio;
		if (verlen < 0.1){
			tmp_ratio = 1.0;
		}
		else{
			tmp_ratio = horlen / verlen;
		} 

		double dist, w;
		double sumW = 0;
		ic::vec3 sumWVal;
		std::map<int, double>::iterator itd;
		bool calc_average = true;
		for (it = nids.begin(); it != nids.end(); ++it) {
			//std::cout << it->first << " ";
			it->second.p.z = it->second.p.z * tmp_ratio * Scale + it->second.p.z * (1.0 - Scale);
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


*/
		double porosity = 1.0;
		if (porType == ic::interpType::CLOUD) {
			// TODO porosity = ic::interpolateScalarTree(PorosityTree, p);
		}
		else if (porType == ic::interpType::SCALAR)
			porosity = porosityValue;
		vel = vel * (1/porosity);

        pp = p;
        vv = vel;

		//if (method == 4){
		//	step = calculate_step(p, vel, lbb, ubb);
		//}

		//auto finish = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> elapsed = finish - start;
		//std::cout << "Calc time Tria: " << elapsed.count() << std::endl;
		//calc_time += elapsed.count();

		count_times++;
		//nids_size = nids.size();
        ic::PrintStat(count_times, FrequencyStat, calc_time, max_calc_time);
	}

/*
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
	*/

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

	void npsatVel::updateStep(double& step){
		step = calculate_step(pp, vv, ll, uu);
	}

	double npsatVel::calculate_step(ic::vec3& p, ic::vec3& v, ic::vec3& l, ic::vec3& u){
		ic::vec3 p1 = p + v * 0.1;
		double t_max = 1000000000;
		double t_min = -1000000000;
		double tmp_min, tmp_max;

		//ic::DEBUG::displayPVasVex(p, v);
		//ic::DEBUG::displayVectorasVex(p1);

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
		//ic::DEBUG::displayVectorasVex(pmin);
		//ic::DEBUG::displayVectorasVex(pmax);
		return pmin.distance(pmax.x, pmax.y, pmax.z)/n_steps;
	}

}

