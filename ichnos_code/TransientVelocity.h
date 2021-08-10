//
// Created by giorg on 8/5/2021.
//

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

namespace TRANS{

    class transVel : public ic::velocityField{
    public:
        transVel(boost::mpi::communicator& world_in, ic::VelType Vtype_in);
        void readVelocityField(std::string vf_file);
        void calcVelocity(ic::vec3& vel, std::map<int, double>& proc_map, ic::vec3& p, double time = 0);
        void reset();
        void updateStep(double& step);

    private:
        bool readXYZfile(std::string prefix, std::string suffix, int ld_zero);
        bool readVelocityFiles(std::string prefix, std::string suffix, int ld_zero);
        bool readVfile(std::string filename, int nPoints, int nSteps, ic::coordDim dim);
        bool readTimeSteps(std::string filename, std::vector<double>& TS);

        ic::VelTR VEL;
        ic::search_tree_trans Tree;

        ic::interpType porType;
        double porosityValue = 1.0;


        int nPoints;
        int nSteps;

        double multiplier = 1.0;
        double Power;
        double Scale = 1.0;
        double initial_diameter = 640;
        double initial_ratio = 20;
        double diameter;
        double ratio;

        double Threshold;
        int FrequencyStat;
        double calc_time = 0.0;
        double max_calc_time = 0.0;
        int count_times = 0;

        double currentTime;

        bool bIsInitialized = false;
        double search_mult = 2.5;
        ic::vec3 ll, uu, pp, vv;
        ic::TimeInterpType timeInterp = ic::TimeInterpType::NEAREST;

    };

    transVel::transVel(boost::mpi::communicator &world_in, ic::VelType Vtype_in)
        :
        velocityField(world_in, Vtype_in)
    {
        InterpolateOutsideDomain = true;
    }

    void transVel::readVelocityField(std::string vf_file) {
        if (world.rank() == 0)
            std::cout << "--> Velocity configuration file: " << vf_file << std::endl;

        po::options_description velocityFieldOptions("Velocity field options");
        po::variables_map vm_vfo;
        velocityFieldOptions.add_options()
            // Velocity parameters
            ("Velocity.Prefix", po::value<std::string>(), "Prefix for the filename")
            ("Velocity.LeadingZeros", po::value<int>()->default_value(4), "e.g 0002->4, 000->3")
            ("Velocity.Suffix", po::value<std::string>(), "ending of file after procid")
            ("Velocity.TimeStepFile", po::value<std::string>(), "This filename with the time steps")
            ("Velocity.TimeInterp", po::value<std::string>(), "Interpolation type between timesteps")
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
        ;

        po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions), vm_vfo);

        { // General
            OwnerThreshold = vm_vfo["General.OwnerThreshold"].as<double>();
            Threshold = vm_vfo["General.Threshold"].as<double>();
            FrequencyStat = vm_vfo["General.FrequencyStat"].as<int>();
        }

        {// Velocity parameters
            multiplier = vm_vfo["Velocity.Multiplier"].as<double>();
            Scale = vm_vfo["Velocity.Scale"].as<double>();
            Power = vm_vfo["Velocity.Power"].as<double>();
            initial_diameter = vm_vfo["Velocity.InitDiameter"].as<double>();
            initial_ratio = vm_vfo["Velocity.InitRatio"].as<double>();

            // Read the time step
            std::string TSfile = vm_vfo["Velocity.TimeStepFile"].as<std::string>();
            std::vector<double> TimeSteps;
            bool tf = readTimeSteps(TSfile, TimeSteps);
            if (!tf) return;

            {
                std::string TimeInterpType = vm_vfo["Velocity.TimeInterp"].as<std::string>();
                if (TimeInterpType.compare("LINEAR") == 0) {
                    timeInterp = ic::TimeInterpType::LINEAR;
                }
            }

            std::string prefix = vm_vfo["Velocity.Prefix"].as<std::string>();
            std::string suffix;
            if (vm_vfo.count("Velocity.Suffix")) {
                suffix = vm_vfo["Velocity.Suffix"].as<std::string>();
            }
            else {
                suffix = ".ich";
            }
            int leadZeros = vm_vfo["Velocity.LeadingZeros"].as<int>();

            // Read the XYZ points
            tf = readXYZfile(prefix, suffix, leadZeros);

            VEL.init(nPoints, nSteps);
            VEL.setTSvalue(TimeSteps);

            // Read the Velocity files
            tf = readVelocityFiles(prefix, suffix, leadZeros);
        }
    }

    bool transVel::readVelocityFiles(std::string prefix, std::string suffix, int ld_zero){
        auto start = std::chrono::high_resolution_clock::now();
        std::string fileVX = prefix + "VX_" + ic::num2Padstr(/*dbg_rank*/world.rank(), ld_zero) + suffix;
        bool tf = readVfile(fileVX, nPoints, nSteps, ic::coordDim::vx);
        if (tf){
            std::string fileVY = prefix + "VY_" + ic::num2Padstr(/*dbg_rank*/world.rank(), ld_zero) + suffix;
            tf = readVfile(fileVY, nPoints, nSteps, ic::coordDim::vy);
        }
        if (tf){
            std::string fileVZ = prefix + "VZ_" + ic::num2Padstr(/*dbg_rank*/world.rank(), ld_zero) + suffix;
            tf = readVfile(fileVZ, nPoints, nSteps, ic::coordDim::vz);
        }
        if (!tf)
            return false;

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read Velocity in : " << elapsed.count() << std::endl;
        return true;
    }

    bool transVel::readVfile(std::string filename, int nPoints, int nSteps, ic::coordDim dim){
        std::cout << "Reading file " + filename << std::endl;
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()){
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }
        else{
            std::string line;
            double v;
            for (int i = 0; i < nPoints; i++){
                getline(datafile, line);
                std::istringstream inp(line.c_str());
                for (int j = 0; j < nSteps; j++){
                    inp >> v;
                    VEL.setVELvalue(v*multiplier, i, j, dim);
                }
            }
        }
        return true;
    }

    bool transVel::readTimeSteps(std::string filename, std::vector<double>& TS){
        std::cout << "Reading file " + filename << std::endl;
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()){
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }
        else{
            std::string line;
            double v;
            while (getline(datafile, line)){
                if (line.size() > 0){
                    std::istringstream inp(line.c_str());
                    inp >> v;
                    TS.push_back(v);
                }
            }
            nSteps = TS.size();
            std::cout << "Proc["<< world.rank() <<  "] Number of time steps: " << nSteps << std::endl;
            datafile.close();
        }
        return true;
    }

    bool transVel::readXYZfile(std::string prefix, std::string suffix, int ld_zero){
        std::string fileXYZ = prefix + "XYZ_" + ic::num2Padstr(/*dbg_rank*/world.rank(), ld_zero) + suffix;
        std::cout << "Reading file " + fileXYZ << std::endl;
        std::ifstream datafile(fileXYZ.c_str());
        if (!datafile.good()){
            std::cout << "Can't open the file " << fileXYZ << std::endl;
            return false;
        }
        else{
            std::vector<ic::cgal_point_3> pp;
            std::vector<ic::TRANS_data> dd;

            std::string line;
            double x, y, z;
            ic::TRANS_data td;
            int cnt = 0;
            while (getline(datafile, line)){
                if (line.size() > 1){
                    std::istringstream inp(line.c_str());
                    inp >> x;
                    inp >> y;
                    inp >> z;
                    inp >> td.proc;
                    inp >> td.diameter;
                    inp >> td.ratio;
                    td.id = cnt;
                    cnt++;

                    pp.push_back(ic::cgal_point_3(x, y, z));
                    dd.push_back(td);
                }
            }
            nPoints = pp.size();
            std::cout << "Proc["<< world.rank() <<  "] Number of XYZ points: " << nPoints << std::endl;
            datafile.close();

            { //Build tree
                auto start = std::chrono::high_resolution_clock::now();
                Tree.insert(boost::make_zip_iterator(boost::make_tuple( pp.begin(),dd.begin() )),
                            boost::make_zip_iterator(boost::make_tuple( pp.end(),dd.end() ) )  );
                Tree.build();
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = finish - start;
                std::cout << "Point Set Building time: " << elapsed.count() << std::endl;
            }
        }
        return true;
    }

    void transVel::calcVelocity(ic::vec3& vel, std::map<int, double>& proc_map, ic::vec3& p, double tm){
        // If this is the first point of this streamline we will carry out one additional range search
        ll.zero();
        uu.zero();
        pp.zero();
        vv.zero();
        if (!bIsInitialized){
            ic::calculate_search_box(p,ll,uu,diameter,ratio,search_mult);
            ic::cgal_point_3 llp(ll.x, ll.y, ll.z);
            ic::cgal_point_3 uup(uu.x, uu.y, uu.z);
            ic::Fuzzy_iso_box_trans fib(llp,uup, 0.0);
            std::vector<boost::tuples::tuple<ic::cgal_point_3, ic::TRANS_data>> tmp;
            Tree.search(std::back_inserter(tmp), fib);
            if (tmp.size() == 0){
                vel = ic::vec3(-99999,-99999,-99999);
                return;
            }
            // Find the closest point
            double mindist = 1000000000;
            double tmp_diam;
            double tmp_ratio;
            ic::TRANS_data closest_point_data;
            for (unsigned int i = 0; i < tmp.size(); ++i) {
                double dist = p.distance(tmp[i].get<0>().x(), tmp[i].get<0>().y(), tmp[i].get<0>().z());
                if (dist < mindist) {
                    mindist = dist;
                    tmp_diam = tmp[i].get<1>().diameter;
                    tmp_ratio = tmp[i].get<1>().ratio;
                }
            }
            diameter = tmp_diam;
            ratio = tmp_ratio;
            bIsInitialized = true;
        }

        // Find the time step
        int i1, i2;
        double t;
        VEL.findIIT(tm, i1, i2, t);

        {
            auto start = std::chrono::high_resolution_clock::now();
            std::map<int, double>::iterator itd;
            std::vector<boost::tuples::tuple<ic::cgal_point_3, ic::TRANS_data>> tmp;
            while (true){
                tmp.clear();
                ic::calculate_search_box(p,ll,uu, diameter, ratio,search_mult);
                ic::cgal_point_3 llp(ll.x, ll.y, ll.z);
                ic::cgal_point_3 uup(uu.x, uu.y, uu.z);
                ic::Fuzzy_iso_box_trans fib(llp,uup, 0.0);
                Tree.search(std::back_inserter(tmp), fib);
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
                    vel = VEL.getVelocity(tmp[i].get<1>().id, i1, i2, t, timeInterp);
                    calc_average = false;
                }
                else{
                    w = 1 / std::pow(scaled_dist, Power);
                    itd = proc_map.find(tmp[i].get<1>().proc);
                    if (itd == proc_map.end()) {
                        proc_map.insert(std::pair<int, double>(tmp[i].get<1>().proc, w));
                    }
                    else{
                        itd->second += w;
                    }
                    sumW += w;
                    sumWVal = sumWVal + VEL.getVelocity(tmp[i].get<1>().id, i1, i2, t, timeInterp) * w;
                }
            }
            if (tmp.size() > 50){
                diameter = tmp_diam;
                ratio = tmp_ratio;
            }

            itd = proc_map.begin();
            for (; itd != proc_map.end(); ++itd) {
                itd->second = itd->second / sumW;
            }

            if (calc_average)
                vel = sumWVal * (1 / sumW);

            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            calc_time += elapsed.count();
            if (elapsed.count() > max_calc_time)
                max_calc_time = elapsed.count();
        }

        double porosity = 1.0;
        if (porType == ic::interpType::CLOUD) {
            // TODO porosity = ic::interpolateScalarTree(PorosityTree, p);
        }
        else if (porType == ic::interpType::SCALAR)
            porosity = porosityValue;
        vel = vel * (1/porosity);

        // TODO Understand why this is needed
        pp = p; // I dont understand this. It seems we need to use these variables in the updateStep method
        vv = vel;

        count_times++;
        ic::PrintStat(count_times, FrequencyStat, calc_time, max_calc_time);

    }

    void transVel::reset() {
        bIsInitialized = false;
        diameter = initial_diameter;
        ratio = initial_ratio;
    }

    void transVel::updateStep(double &step) {

    }
}

