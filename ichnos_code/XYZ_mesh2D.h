//
// Created by giorg on 11/5/2021.
//

#ifndef ICHNOS_XYZ_MESH2D_H
#define ICHNOS_XYZ_MESH2D_H

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include "ichnos_structures.h"
#include "ichnos_XYZ_base.h"
#include "ichnos_mesh.h"

namespace po = boost::program_options;

namespace ICHNOS{
    class XYZ_MESH2D : public XYZ_base {
    public:
        XYZ_MESH2D(boost::mpi::communicator& world_in);
        bool readXYZdata(std::string vf_file);
        void calcWeights(vec3& p,
                         std::vector<int>& ids,
                         std::vector<double>& weights,
                         std::map<int, double>& proc_map,
                         helpVars& pvlu, bool& out);
        void reset(Streamline& S);
        int getINTInfo(infoType I);
        void sendVec3Data(std::vector<vec3>& data){};
        int getNpnts(){return static_cast<int>(mesh.size());}
        void SetInterpType(MeshVelInterpType tp){velInterpType = tp;};

    private:
        int nLay;
        int nElements;
        //double diameter;
        double initial_diameter = 0.0;

        //Mesh2D MSH;
        std::vector<std::vector<double>> elev;
        std::vector<std::vector<int>> mesh;
        std::vector<vec3> nodes;
        search_tree_iwfm Tree;

        bool readNodes(std::string filename);
        bool readElev(std::string filename);
        bool readMesh(std::string filename);
        //bool readFaceIds(std::string filename);
        MeshVelInterpType velInterpType = MeshVelInterpType::UNKNOWN;

        const double Threshold = 0.01;
        const double Power = 3;
        const double invTranfTol = 0.0001;
    };

    XYZ_MESH2D::XYZ_MESH2D(boost::mpi::communicator& world_in)
            :
            XYZ_base(world_in)
    {}

    bool XYZ_MESH2D::readXYZdata(std::string vf_file){
        if (world.rank() == 0){
            std::cout << "\tReading XYZ data from " << vf_file << "..." << std::endl;
        }

        po::options_description velocityFieldOptions("Velocity field options");
        po::variables_map vm_vfo;
        velocityFieldOptions.add_options()
            ("MESH2D.NodeFile", po::value<std::string>(), "An array of the node coordinates")
            ("MESH2D.Meshfile", po::value<std::string>(), "An array of the Mesh2D ids")
            ("MESH2D.ElevationFile", po::value<std::string>(), "An array of the Elevations")
            ("MESH2D.Nlayers", po::value<int>()->default_value(4), "Number of layers")
            ("Velocity.LeadingZeros", po::value<int>()->default_value(4), "e.g 0002->4, 000->3")
            ("Velocity.Suffix", po::value<std::string>(), "ending of file after proc id")

            //("MESH2D.SearchDiameter", po::value<double>()->default_value(5000), "The search diameter")
        //("MESH2D.FaceIdFile", po::value<std::string>(), "Face ids for each element")
            ;
        po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions,true), vm_vfo);

        std::string nodefile = vm_vfo["MESH2D.NodeFile"].as<std::string>();
        std::string Meshfile = vm_vfo["MESH2D.Meshfile"].as<std::string>();
        std::string Elevfile = vm_vfo["MESH2D.ElevationFile"].as<std::string>();

        nLay = vm_vfo["MESH2D.Nlayers"].as<int>();
        //initial_diameter = vm_vfo["MESH2D.SearchDiameter"].as<double>();
        //initial_diameter = initial_diameter*initial_diameter;

        if (world.size() > 1 && runAsThread == 0){
            int leadingZeros = vm_vfo["Velocity.LeadingZeros"].as<int>();
            std::string suffix = vm_vfo["Velocity.Suffix"].as<std::string>();
            nodefile = nodefile + ic::num2Padstr(world.rank(), leadingZeros) + suffix;
            Elevfile = Elevfile + ic::num2Padstr(world.rank(), leadingZeros) + suffix;
            Meshfile = Meshfile + ic::num2Padstr(world.rank(), leadingZeros) + suffix;
        }

        bool tf = readNodes(nodefile);
        if (!tf){return false;}
        tf = readElev(Elevfile);
        if (!tf){return false;}
        tf = readMesh(Meshfile);

        if (world.size() > 1 && runAsThread == 0){
            std::cout << world.rank() << "ND: " << nodes.size() << ", MSH: " << mesh.size() << ", ELEV: " << elev.size() << std::endl;
        }

        return true;
    }

    bool XYZ_MESH2D::readNodes(std::string filename) {
        std::vector<std::vector<double>> data;
        bool tf = READ::read2Darray<double>(filename, 2, data);
        if (!tf){return false;}
        for (unsigned int i = 0; i < data.size(); ++i){
            nodes.push_back(vec3(data[i][0], data[i][1], 0.0));
        }
        //MSH.setNodes(nodes);
        return true;
    }

    bool XYZ_MESH2D::readElev(std::string filename) {
        return READ::read2Darray<double>(filename, nLay+1, elev);
    }

    /*
    bool XYZ_MESH2D::readFaceIds(std::string filename) {
        std::vector<std::vector<int>> data;
        bool tf = READ::read2Darray<int>(filename, 4, data);
        if (!tf){return false;}

        for (unsigned int i = 0; i < data.size(); ++i){
            int nvert = 4;
            if (data[i][3] == 0)
                nvert = 3;
            std::vector<int> ids;
            for (int ii = 0; ii < nvert; ++ii){
                int id = boost::math::sign(data[i][ii]) * (abs(data[i][ii]) - 1);
                ids.push_back(id);
            }
            faceIds.push_back(ids);
        }
        nFaces = faceIds.size();
        return true;
    }
    */

    bool XYZ_MESH2D::readMesh(std::string filename){
        std::vector<std::vector<int>> data;
        bool tf = READ::read2Darray<int>(filename, 5, data);
        if (!tf){return false;}

        Pnt_IWFM_info td;
        std::vector<cgal_point_3> pp;
        std::vector<Pnt_IWFM_info> dd;

        initial_diameter = 0.0;
        for (unsigned int i = 0; i < data.size(); ++i){
            int nvert = 4;
            if (data[i][3] == 0)
                nvert = 3;
            double diam;
            vec3 elem_bc;
            { // Calculate element barycenter and diameter
                for (int ii = 0; ii < nvert; ++ii){
                    elem_bc.x += nodes[data[i][ii]-1].x;
                    elem_bc.y += nodes[data[i][ii]-1].y;
                    elem_bc.z += nodes[data[i][ii]-1].z;
                }
                double n = static_cast<double>(nvert);
                elem_bc.x = elem_bc.x / n;
                elem_bc.y = elem_bc.y / n;
                elem_bc.z = elem_bc.z / n;

                // Calculate average distance between element barycenter and element nodes
                // Set the Initial search diameter equal to 2 times the largest diameter
                double dst = 0;

                for (int ii = 0; ii < nvert; ++ii){
                    dst += elem_bc.distance(nodes[data[i][ii]-1].x, nodes[data[i][ii]-1].y, nodes[data[i][ii]-1].z);
                }
                diam = 2.0*(dst/n);
                if (diam > initial_diameter){
                    initial_diameter = diam;
                }

            }

            std::vector<int> ids;
            for (int ii =0; ii < nvert; ++ii){
                ids.push_back(data[i][ii]-1);
            }
            pp.push_back(cgal_point_3(elem_bc.x, elem_bc.y, elem_bc.z));
            td.elem_id = i;
            td.proc = data[i][4];
            td.tri_id = mesh.size();
            td.diameter = diam;
            mesh.push_back(ids);

            dd.push_back(td);
        }
        initial_diameter = initial_diameter*2.0;
        nElements = mesh.size();

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

    void XYZ_MESH2D::calcWeights(vec3 &p,
                                 std::vector<int> &ids,
                                 std::vector<double> &weights,
                                 std::map<int, double>& proc_map,
                                 helpVars& pvlu, bool& out){
        out = false;
        ids.clear();
        weights.clear();
        vec3 p2D (p.x, p.y, 0.0);

        std::map<int, double>::iterator itd;
        std::vector<boost::tuples::tuple<cgal_point_3, Pnt_IWFM_info>> tmp;
        while (true){
            tmp.clear();
            Fuzzy_sphere_iwfm fc(cgal_point_3(p2D.x, p2D.y, p2D.z), pvlu.diameter);
            Tree.search(std::back_inserter(tmp), fc);
            if (tmp.size() >= 1){
                break;
            }
            else{
                pvlu.diameter = pvlu.diameter*1.5;
                if (pvlu.diameter > initial_diameter){
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
        vec3 bc; // iso parametric coordinates
        int iLay;
        int elem_id;
        double w;
        double sumW = 0;
        for (unsigned int i = 0; i < tmp.size(); ++i){
            if (!elem_found){
                int el_id = tmp[i].get<1>().tri_id;
                if (mesh[el_id].size() == 3){
                    bool tf = isPointInTriangle(p2D,nodes[mesh[el_id][0]],nodes[mesh[el_id][1]],nodes[mesh[el_id][2]], bc);
                    if (tf){
                        elem_found = true;
                    }
                }
                else{
                    bool tf = isPointInTriangle(p2D,nodes[mesh[el_id][0]],nodes[mesh[el_id][1]],nodes[mesh[el_id][2]], bc);
                    if (!tf){
                        tf = isPointInTriangle(p2D,nodes[mesh[el_id][0]],nodes[mesh[el_id][2]],nodes[mesh[el_id][3]], bc);
                    }
                    if (tf){
                        QuadInverseMappingV1(p2D, nodes[mesh[el_id][0]], nodes[mesh[el_id][1]],
                                             nodes[mesh[el_id][2]], nodes[mesh[el_id][3]], bc);
                        //QuadInverseMapping(p2D, nodes[mesh[el_id][0]], nodes[mesh[el_id][1]],
                        //                            nodes[mesh[el_id][2]], nodes[mesh[el_id][3]], bc, invTranfTol);
                        if (tf){
                            elem_found = true;
                        }
                    }
                }
                if (elem_found){
                    elem_id = el_id;
                }
            }

            // Calculate the distance between the particle and the tree points
            double dist = p2D.distance(tmp[i].get<0>().x(), tmp[i].get<0>().y(), tmp[i].get<0>().z());
            if (dist < Threshold){
                proc_map.clear();
                proc_map.insert(std::pair<int, double>(tmp[i].get<1>().proc, 1.0));
                if (elem_found)
                    break;
            }
            else{
                w = 1 / std::pow(dist, Power);
                itd = proc_map.find(tmp[i].get<1>().proc);
                if (itd == proc_map.end()){
                    proc_map.insert(std::pair<int, double>(tmp[i].get<1>().proc, w));
                }
                else{
                    itd->second += w;
                }
                sumW += w;
            }
        }

        if (elem_found){
            // Find the layer
            double ztop;
            double zbot;
            for (int i = 0; i < nLay; ++i){
                if (mesh[elem_id].size() == 3){
                    ztop = elev[mesh[elem_id][0]][i]*bc.x + elev[mesh[elem_id][1]][i]*bc.y +elev[mesh[elem_id][2]][i]*(1-bc.x-bc.y);
                    zbot = elev[mesh[elem_id][0]][i+1]*bc.x + elev[mesh[elem_id][1]][i+1]*bc.y +elev[mesh[elem_id][2]][i+1]*(1-bc.x-bc.y);
                }
                else if (mesh[elem_id].size() == 4){
                    double n1, n2, n3, n4;
                    QuadShapeFunctions(bc.x, bc.y, n1, n2, n3, n4);
                    ztop = elev[mesh[elem_id][0]][i]*n1 + elev[mesh[elem_id][1]][i]*n2 + elev[mesh[elem_id][2]][i]*n3 + elev[mesh[elem_id][3]][i]*n4;
                    zbot = elev[mesh[elem_id][0]][i+1]*n1 + elev[mesh[elem_id][1]][i+1]*n2 + elev[mesh[elem_id][2]][i+1]*n3 + elev[mesh[elem_id][3]][i+1]*n4;
                }
                if (i == 0){
                    if (p.z >= zbot){
                        iLay = i;
                        out = true;
                        break;
                    }
                }
                else if (i == nLay - 1){
                    if (p.z < ztop){
                        iLay = i;
                        out = true;
                        break;
                    }
                }
                else{
                    if (p.z >= zbot && p.z < ztop){
                        iLay = i;
                        out = true;
                        break;
                    }
                }
            }
            // parametric value of the z coodrinate
            bc.z = (p.z - zbot)/(ztop - zbot);

            if (velInterpType == MeshVelInterpType::ELEMENT){
                weights.push_back(1.0);
            }
            else{
                weights.push_back(bc.x);
                weights.push_back(bc.y);
                weights.push_back(bc.z);
            }
            ids.push_back(elem_id);
            ids.push_back(iLay);
            itd = proc_map.begin();
            for (; itd != proc_map.end(); ++itd) {
                itd->second = itd->second / sumW;
            }

            pvlu.ll = 999999999999;
            pvlu.uu = -999999999999;
            for (unsigned int i = 0; i < mesh[elem_id].size(); ++i){
                if (pvlu.ll.x > nodes[mesh[elem_id][i]].x){pvlu.ll.x = nodes[mesh[elem_id][i]].x;}
                if (pvlu.ll.y > nodes[mesh[elem_id][i]].y){pvlu.ll.y = nodes[mesh[elem_id][i]].y;}
                if (pvlu.ll.z > elev[mesh[elem_id][i]][iLay+1]){pvlu.ll.z = elev[mesh[elem_id][i]][iLay+1];}
                if (pvlu.uu.x < nodes[mesh[elem_id][i]].x){pvlu.uu.x = nodes[mesh[elem_id][i]].x;}
                if (pvlu.uu.y < nodes[mesh[elem_id][i]].y){pvlu.uu.y = nodes[mesh[elem_id][i]].y;}
                if (pvlu.uu.z < elev[mesh[elem_id][i]][iLay]){pvlu.uu.z = elev[mesh[elem_id][i]][iLay];}
            }
        }
        else{
            weights.push_back(1.0);
            ids.push_back(-9);
            ids.push_back(-9);
            pvlu.ll.zero();
            pvlu.uu.zero();
        }
    }

    void XYZ_MESH2D::reset(Streamline& S) {
        S.PVLU.diameter = initial_diameter;
    }
    int XYZ_MESH2D::getINTInfo(ICHNOS::infoType I) {
        if (I == infoType::Nelem){
            return nElements;
        }
        return 0;
    }

}

#endif //ICHNOS_XYZ_MESH2D_H
