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
        virtual bool readXYZdata(std::string vf_file){return true;}
        virtual void calcWeights(vec3& p,
                                 std::vector<int>& ids,
                                 std::vector<double>& weights,
                                 std::map<int, double>& proc_map,
                                 bool& out){}
        virtual void reset(){}
        virtual void sendVec3Data(std::vector<vec3>& data){};

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
        bool readXYZdata(std::string vf_file);
        void calcWeights(vec3& p,
                         std::vector<int>& ids,
                         std::vector<double>& weights,
                         std::map<int, double>& proc_map,
                         bool& out);
        void reset();
        void sendVec3Data(std::vector<vec3>& data);
        int getNpnts(){return Tree.size();}
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

    bool XYZ_cloud::readXYZdata(std::string vf_file) {
        if (world.rank() == 0){
            std::cout << "\tReading XYZ data..." << std::endl;
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
            return false;
        }
        std::string fileXYZ;

        if (vtype == VelType::STEADY || suffix.compare(".h5") == 0){
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
        if (!tf)
            return false;

        { //Build tree
            auto start = std::chrono::high_resolution_clock::now();
            Tree.insert(boost::make_zip_iterator(boost::make_tuple( pp.begin(),dd.begin() )),
                        boost::make_zip_iterator(boost::make_tuple( pp.end(),dd.end() ) )  );
            Tree.build();
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            std::cout << "\tPoint Set Building time: " << elapsed.count() << std::endl;
        }
        return true;
    }

    void XYZ_cloud::calcWeights(vec3 &p,
                                std::vector<int> &ids,
                                std::vector<double> &weights,
                                std::map<int, double>& proc_map,
                                bool& out) {
        ids.clear();
        weights.clear();
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
                        std::cout << "I cant find any point around ("
                                  << p.x << "," << p.y << "," << p.z
                                  <<") within the initial diameter of " << initial_diameter
                                  << ". Consider increasing the Initial diameter" << std::endl;
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
                    proc_map.clear();
                    proc_map.insert(std::pair<int, double>(tmp[i].get<1>().proc, 1.0));
                    break;
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

    void XYZ_cloud::sendVec3Data(std::vector<vec3> &data) {
        data.clear();
        data.push_back(pp);
        data.push_back(ll);
        data.push_back(uu);
    }

    class XYZ_IWFM : public XYZ_base {
    public:
        XYZ_IWFM(boost::mpi::communicator& world_in);
        bool readXYZdata(std::string vf_file);
        void calcWeights(vec3& p,
                         std::vector<int>& ids,
                         std::vector<double>& weights,
                         std::map<int, double>& proc_map,
                         bool& out);
        void reset();
        void sendVec3Data(std::vector<vec3>& data);
        int getNpnts(){return nXYpnts * nLay;}

    private:
        search_tree_iwfm Tree;
        VelType vtype;
        std::vector<std::vector<double>> elev;
        std::vector<std::vector<int>> tria_mesh;
        std::vector<vec3> nodes;
        double diameter;
        double initial_diameter = 640;
        double Threshold;
        double Power;
        int nLay;
        int nXYpnts;
        double search_mult = 1.5;

        bool readNodes(std::string filename);
        bool readStrat(std::string filename);
        bool readMesh(std::string filename);
        vec3 ll, uu, pp, vv;
    };
    XYZ_IWFM::XYZ_IWFM(boost::mpi::communicator& world_in)
        :
        XYZ_base(world_in)
    {}

    bool XYZ_IWFM::readXYZdata(std::string vf_file) {
        if (world.rank() == 0){
            std::cout << "\tReading XYZ data from " << vf_file << "..." << std::endl;
        }

        po::options_description velocityFieldOptions("Velocity field options");
        po::variables_map vm_vfo;
        velocityFieldOptions.add_options()
            ("IWFM.NodeFile", po::value<std::string>(), "An array of the node coordinates")
            ("IWFM.ElevationFile", po::value<std::string>(), "An array of the Elevations")
            ("IWFM.Meshfile", po::value<std::string>(), "An array of the Mesh ids")
            ("IWFM.Nlayers", po::value<int>()->default_value(4), "Number of layers")
            ("IWFM.Power", po::value<double>()->default_value(3.0), "Power of the IDW interpolation")
            ("IWFM.Diameter", po::value<double>()->default_value(1.0), "Scale the domain before velocity calculation")
            ("IWFM.Threshold", po::value<double>()->default_value(5000), "Initial diameter")
            ;

        po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions,true), vm_vfo);

        std::string nodefile = vm_vfo["IWFM.NodeFile"].as<std::string>();
        std::string Stratfile = vm_vfo["IWFM.ElevationFile"].as<std::string>();
        std::string Meshfile = vm_vfo["IWFM.Meshfile"].as<std::string>();
        nLay = vm_vfo["IWFM.Nlayers"].as<int>();
        Threshold = vm_vfo["IWFM.Threshold"].as<double>();
        initial_diameter = vm_vfo["IWFM.Diameter"].as<double>();
        Power = vm_vfo["IWFM.Power"].as<double>();
        diameter = initial_diameter;

        bool tf = readNodes(nodefile);
        if (!tf){return false;}
        tf = readStrat(Stratfile);
        if (!tf){return false;}
        tf = readMesh(Meshfile);

        return true;
    }

    void XYZ_IWFM::calcWeights(vec3 &p,
                               std::vector<int> &ids,
                               std::vector<double> &weights,
                               std::map<int, double>& proc_map,
                               bool& out) {
        out = false;
        ids.clear();
        weights.clear();
        vec3 p2D (p.x, p.y, 0.0);

        std::map<int, double>::iterator itd;
        std::vector<boost::tuples::tuple<cgal_point_3, Pnt_IWFM_info>> tmp;
        while (true){
            tmp.clear();
            Fuzzy_sphere_iwfm fc(cgal_point_3(p2D.x, p2D.y, p2D.z), diameter);
            Tree.search(std::back_inserter(tmp), fc);
            if (tmp.size() >= 1){
                break;
            }
            else{
                diameter = diameter*1.5;
                if (diameter > initial_diameter){
                    std::cout << "I cant find any point around ("
                              << p2D.x << "," << p2D.y << "," << p2D.z
                              <<") within the initial diameter of " << initial_diameter
                              << ". Consider increasing the Initial diameter" << std::endl;
                    out = false;
                    return;
                }
            }
        }
        bool elem_found = false;
        double w;
        double sumW = 0;
        double diam, ratio;
        int elem_id, tri_id, near_id;
        double closest_dist = 1000000000;
        vec3 nearest_bc;
        for (unsigned int i = 0; i < tmp.size(); ++i){
            vec3 bc;
            if (!elem_found){
                bool tf = isPointInTriangle(p2D,
                                            nodes[tria_mesh[tmp[i].get<1>().tri_id][0]],
                                            nodes[tria_mesh[tmp[i].get<1>().tri_id][1]],
                                            nodes[tria_mesh[tmp[i].get<1>().tri_id][2]],
                                            bc);
                if (tf){
                    near_id = i;
                    elem_id = tmp[i].get<1>().elem_id;
                    tri_id = tmp[i].get<1>().tri_id;
                    diam = tmp[i].get<1>().diameter;
                    nearest_bc = bc;
                    elem_found = true;
                }
            }

            // Calculate the distance
            double dist = p2D.distance(tmp[i].get<0>().x(), tmp[i].get<0>().y(), tmp[i].get<0>().z());
            if (dist < closest_dist && !elem_found){
                closest_dist = dist;
                near_id = i;
                elem_id = tmp[i].get<1>().elem_id;
                tri_id = tmp[i].get<1>().tri_id;
                diam = tmp[i].get<1>().diameter;
                nearest_bc = bc;
            }
            if (dist < Threshold){
                proc_map.clear();
                proc_map.insert(std::pair<int, double>(tmp[i].get<1>().proc, 1.0));
                if (elem_found)
                    break;
            }
            else{
                w = 1 / std::pow(dist, Power);
                itd = proc_map.find(tmp[i].get<1>().proc);
                if (itd == proc_map.end()) {
                    proc_map.insert(std::pair<int, double>(tmp[i].get<1>().proc, w));
                }
                else{
                    itd->second += w;
                }
                sumW += w;
            }
        }
        //if (!elem_found){
        //    std::cout << "elem not found" << std::endl;
        //}

        {
            // Find the layer that the particle is in
            double lay_t = elev[tria_mesh[tri_id][0]][1]*nearest_bc.x + elev[tria_mesh[tri_id][1]][1]*nearest_bc.y + elev[tria_mesh[tri_id][2]][1]*nearest_bc.z;
            if (p.z >= lay_t){
                ids.push_back(elem_id);
                ratio = diam/tmp[near_id].get<1>().height[0];
            }
            else{
                for (int j = 2; j < nLay-1; ++j){
                    double lay_b = elev[tria_mesh[tri_id][0]][j]*nearest_bc.x + elev[tria_mesh[tri_id][1]][j]*nearest_bc.y + elev[tria_mesh[tri_id][2]][j]*nearest_bc.z;
                    if (j == nLay-1){
                        if (p.z < lay_b){
                            ids.push_back(elem_id + nXYpnts*(nLay-1));
                            ratio = diam/tmp[near_id].get<1>().height[nLay-1];
                        }
                    }
                    else{
                        if (p.z <= lay_t && p.z >= lay_b){
                            ids.push_back(elem_id + (j-1)*nXYpnts);
                            ratio = diam/tmp[near_id].get<1>().height[j-1];
                            break;
                        }
                        else{
                            lay_t = lay_b;
                        }
                    }
                }
            }
            weights.push_back(1.0);
            elem_found = true;
            out = true;
        }


        itd = proc_map.begin();
        for (; itd != proc_map.end(); ++itd) {
            itd->second = itd->second / sumW;
        }
        calculate_search_box(p, ll, uu, diam, ratio, search_mult);
        pp = p;
    }

    void XYZ_IWFM::reset() {

    }

    void XYZ_IWFM::sendVec3Data(std::vector<vec3> &data) {
        data.clear();
        data.push_back(pp);
        data.push_back(ll);
        data.push_back(uu);
    }

    bool XYZ_IWFM::readNodes(std::string filename) {
        std::vector<std::vector<double>> data;
        bool tf = READ::read2Darray<double>(filename, 2, data);
        if (!tf){return false;}
        for (unsigned int i = 0; i < data.size(); ++i){
            nodes.push_back(vec3(data[i][0], data[i][1], 0.0));
        }
        return true;
    }

    bool XYZ_IWFM::readStrat(std::string filename) {
        return READ::read2Darray<double>(filename, nLay+1, elev);
    }

    bool XYZ_IWFM::readMesh(std::string filename) {
        std::vector<std::vector<int>> data;
        bool tf = READ::read2Darray<int>(filename, 5, data);
        if (!tf){return false;}
        int i1, i2, i3;
        double oneThird = 1.0/3.0;
        double oneForth = 1.0/4.0;
        Pnt_IWFM_info td;
        std::vector<cgal_point_3> pp;
        std::vector<Pnt_IWFM_info> dd;
        nXYpnts = data.size();
        for (unsigned int i = 0; i < data.size(); ++i){
            int nvert = 4;
            if (data[i][3] == 0)
                nvert = 3;
            double diam;
            std::vector<double> av_heights;
            {// Calculate diameter and element height
                // Calculate element barycenter
                vec3 elem_bc;
                for (int ii = 0; ii < nvert; ++ii){
                    elem_bc = elem_bc + nodes[data[i][ii]-1];
                }
                elem_bc = elem_bc * (1/static_cast<double>(nvert));

                // Calculate average distance between element barycenter and element nodes
                double dst = 0;
                for (int ii = 0; ii < nvert; ++ii){
                    dst += elem_bc.distance(nodes[data[i][ii]-1].x, nodes[data[i][ii]-1].y, nodes[data[i][ii]-1].z);
                }
                diam = 2.0*(dst/static_cast<double>(nvert));

                // Calculate average layer heights
                for (int j = 0; j < nLay; ++j){
                    double hght = 0;
                    for (int k = 0; k < nvert; ++k){
                        hght += elev[data[i][k]-1][j] - elev[data[i][k]-1][j+1];
                    }
                    av_heights.push_back(hght/static_cast<double>(nvert));
                }
            }

            if (data[i][3] == 0){
                // The element is a triangle
                std::vector<int> ids;
                i1 = data[i][0]-1; ids.push_back(i1);
                i2 = data[i][1]-1; ids.push_back(i2);
                i3 = data[i][2]-1; ids.push_back(i3);
                vec3 bc = (nodes[i1] + nodes[i2] + nodes[i3])*oneThird;
                pp.push_back(cgal_point_3(bc.x, bc.y, bc.z));
                td.elem_id = i;
                td.proc = data[i][4];
                td.tri_id = tria_mesh.size();
                td.diameter = diam;
                td.height = av_heights;

                tria_mesh.push_back(ids);
                dd.push_back(td);

            }
            else{
                // The element is a quadrilateral. We have to split it to 2 triangles
                // 1st triangle
                std::vector<int> ids;
                i1 = data[i][0]-1; ids.push_back(i1);
                i2 = data[i][1]-1; ids.push_back(i2);
                i3 = data[i][2]-1; ids.push_back(i3);
                vec3 bc = (nodes[i1] + nodes[i2] + nodes[i3])*oneThird;
                pp.push_back(cgal_point_3(bc.x, bc.y, bc.z));
                td.elem_id = i;
                td.proc = data[i][4];
                td.tri_id = tria_mesh.size();
                td.diameter = diam;
                td.height = av_heights;
                dd.push_back(td);
                tria_mesh.push_back(ids);

                // 2nd triangle
                std::vector<int> ids1;
                i1 = data[i][0]-1; ids1.push_back(i1);
                i2 = data[i][2]-1; ids1.push_back(i2);
                i3 = data[i][3]-1; ids1.push_back(i3);
                bc = (nodes[i1] + nodes[i2] + nodes[i3])*oneThird;
                pp.push_back(cgal_point_3(bc.x, bc.y, bc.z));
                td.elem_id = static_cast<int>(i);
                td.proc = data[i][4];
                td.tri_id = tria_mesh.size();
                td.diameter = diam;
                td.height = av_heights;
                dd.push_back(td);
                tria_mesh.push_back(ids1);
            }
        }

        {//Build tree
            auto start = std::chrono::high_resolution_clock::now();
            Tree.insert(boost::make_zip_iterator(boost::make_tuple( pp.begin(),dd.begin() )),
                        boost::make_zip_iterator(boost::make_tuple( pp.end(),dd.end() ) )  );
            Tree.build();
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            std::cout << "\tPoint Set Building time: " << elapsed.count() << std::endl;
        }
        return true;
    }
}
