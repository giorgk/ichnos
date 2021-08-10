#pragma once
#include <chrono>
#include <ctime>
#include <utility>

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include "ichnos_structures.h"

namespace po = boost::program_options;

namespace ICHNOS{

    class XYZ_base{
    public:
        XYZ_base(boost::mpi::communicator& world_in);
        virtual void readXYZdata(std::string vf_file){}
        virtual void calcWeights(vec3& p,
                                 std::vector<int>& ids,
                                 std::vector<double>& weights,
                                 std::map<int, double>& proc_map,
                                 bool& out){}
        virtual void reset(){}

    protected:
        boost::mpi::communicator world;
    };

    XYZ_base::XYZ_base(boost::mpi::communicator& world_in)
        :
        world(world_in)
    {}


    class XYZ_cloud : public XYZ_base {
    public:
        XYZ_cloud(boost::mpi::communicator& world_in);
        void readXYZdata(std::string vf_file);
        void calcWeights(vec3& p,
                         std::vector<int>& ids,
                         std::vector<double>& weights,
                         std::map<int, double>& proc_map,
                         bool& out);
        void reset();
    private:
        search_tree_info Tree;
        VelType vtype;

        double Power;
        double Scale = 1.0;
        double diameter;
        double ratio;
        double Threshold;
        double initial_diameter = 640;
        double initial_ratio = 20;

        bool bIsInitialized = false;
        double search_mult = 2.5;
        vec3 ll, uu, pp, vv;
    };

    XYZ_cloud::XYZ_cloud(boost::mpi::communicator& world_in)
        :
        XYZ_base(world_in)
    {}

    void XYZ_cloud::readXYZdata(std::string vf_file) {
        if (world.rank() == 0){
            std::cout << "Reading XYZ data..." << std::endl;
        }

        po::options_description velocityFieldOptions("Velocity field options");
        po::variables_map vm_vfo;
        velocityFieldOptions.add_options()
            ("Velocity.Prefix", po::value<std::string>(), "Prefix for the filename")
            ("Velocity.LeadingZeros", po::value<int>()->default_value(4), "e.g 0002->4, 000->3")
            ("Velocity.Suffix", po::value<std::string>(), "ending of file after procid")
            ("Velocity.Type", po::value<std::string>(), "Type of velocity. (STEADY or TRANS)")
            ("Velocity.Power", po::value<double>()->default_value(3.0), "Power of the IDW interpolation")
            ("Velocity.Scale", po::value<double>()->default_value(1.0), "Scale the domain before velocity calculation")
            ("Velocity.InitDiameter", po::value<double>()->default_value(5000), "Initial diameter")
            ("Velocity.InitRatio", po::value<double>()->default_value(1), "Initial ratio")
            ("General.Threshold", po::value<double>()->default_value(0.001), "Threshold of distance of IDW")
            ;

        po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions,true), vm_vfo);

        std::string prefix = vm_vfo["Velocity.Prefix"].as<std::string>();
        std::string suffix;
        if (vm_vfo.count("Velocity.Suffix")) {
            suffix = vm_vfo["Velocity.Suffix"].as<std::string>();
        }
        else {
            suffix = ".ich";
        }
        int leadZeros = vm_vfo["Velocity.LeadingZeros"].as<int>();
        vtype = castVelType2Enum(vm_vfo["Velocity.Type"].as<std::string>());
        if (vtype == VelType::INVALID) {
            std::cout << vm_vfo["Velocity.Type"].as<std::string>() << " is an invalid velocity type" << std::endl;
            return;
        }
        std::string fileXYZ;
        if (vtype == VelType::STEADY){
            fileXYZ = prefix + num2Padstr(/*dbg_rank*/world.rank(), leadZeros) + suffix;
        }
        else if (vtype == VelType::TRANS){
            fileXYZ = prefix + "XYZ_" + num2Padstr(/*dbg_rank*/world.rank(), leadZeros) + suffix;
        }

        Threshold = vm_vfo["General.Threshold"].as<double>();
        Scale = vm_vfo["Velocity.Scale"].as<double>();
        Power = vm_vfo["Velocity.Power"].as<double>();
        initial_diameter = vm_vfo["Velocity.InitDiameter"].as<double>();
        initial_ratio = vm_vfo["Velocity.InitRatio"].as<double>();

        std::vector<cgal_point_3> pp;
        std::vector<Pnt_info> dd;
        bool tf = READ::readXYZfile(fileXYZ, pp, dd);

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

    void XYZ_cloud::calcWeights(vec3 &p,
                                std::vector<int> &ids,
                                std::vector<double> &weights,
                                std::map<int, double>& proc_map,
                                bool& out) {
        ll.zero();
        uu.zero();
        pp.zero();
        vv.zero();
        if (!bIsInitialized){
            calculate_search_box(p,ll,uu,diameter,ratio,search_mult);
            cgal_point_3 llp(ll.x, ll.y, ll.z);
            cgal_point_3 uup(uu.x, uu.y, uu.z);
            Fuzzy_iso_box_info fib(llp, uup, 0.0);
            std::vector<boost::tuples::tuple<cgal_point_3, Pnt_info>> tmp;
            Tree.search(std::back_inserter(tmp), fib);
            if (tmp.size() == 0){
                out = false;
                return;
            }
            // Find the closest point
            double mindist = 1000000000;
            double tmp_diam;
            double tmp_ratio;
            Pnt_info closest_point_data;
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

        {
            std::map<int, double>::iterator itd;
            std::vector<boost::tuples::tuple<cgal_point_3, Pnt_info>> tmp;
            while (true){
                tmp.clear();
                calculate_search_box(p,ll,uu, diameter, ratio,search_mult);
                cgal_point_3 llp(ll.x, ll.y, ll.z);
                cgal_point_3 uup(uu.x, uu.y, uu.z);
                Fuzzy_iso_box_info fib(llp, uup, 0.0);
                Tree.search(std::back_inserter(tmp), fib);
                if (tmp.size() >= 3){
                    break;
                }
                else{
                    diameter = diameter*1.5;
                    if (diameter > initial_diameter){
                        out = false;
                        return;
                    }
                }
            }

            // por_xyz refers to point translated to origin
            double porx, pory, porz, scaled_dist, actual_dist, w;
            bool calc_average = true;
            double sumW = 0;
            vec3 sumWVal;
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
                    ids.clear();
                    weights.clear();
                    ids.push_back(tmp[i].get<1>().id);
                    weights.push_back(1.0);
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
                    ids.push_back(tmp[i].get<1>().id);
                    weights.push_back(w);
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
            pp = p;
        }
    }
    void XYZ_cloud::reset() {
        bIsInitialized = false;
        diameter = initial_diameter;
        ratio = initial_ratio;
    }

    class XYZ_IWFM : public XYZ_base {
    public:
        XYZ_IWFM(boost::mpi::communicator& world_in);
        void readXYZdata(std::string vf_file);
        void calcWeights(vec3& p,
                         std::vector<int>& ids,
                         std::vector<double>& weights,
                         std::map<int, double>& proc_map,
                         bool& out);
        void reset();
    };
    XYZ_IWFM::XYZ_IWFM(boost::mpi::communicator& world_in)
        :
        XYZ_base(world_in)
    {}

    void XYZ_IWFM::readXYZdata(std::string vf_file) {
        std::cout << "IWFM version" << std::endl;
        std::cout << "File: " << vf_file << std::endl;
    }

    void XYZ_IWFM::calcWeights(vec3 &p,
                               std::vector<int> &ids,
                               std::vector<double> &weights,
                               std::map<int, double>& proc_map,
                               bool& out) {
        std::cout << "IWFM version" << std::endl;
    }

    void XYZ_IWFM::reset() {
        std::cout << "CLOUD version" << std::endl;
    }
}
