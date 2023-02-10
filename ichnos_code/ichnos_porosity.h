//
// Created by giorgk on 2/9/2023.
//

#ifndef ICHNOS_ICHNOS_POROSITY_H
#define ICHNOS_ICHNOS_POROSITY_H

#include "ichnos_structures.h"
#include "ichnos_utils.h"


namespace ICHNOS{
    class Porosity_base {
    public:
        Porosity_base(){};
        bool readData(std::string userInput);
        double interpolate(vec3 p);
        void reset();

    private:
        interpType porInterpolation;
        double porValue = 1.0;
        bool readFile(std::string filename);
        search_tree_info Tree;
        std::vector<double> PorosityArray;
        double cloudInterp(vec3 p);

        double searchDiameter;
        double ratioCurrent;
        double ratioInit;

        double Radius;
        double Power;
        double Threshold;
    };

    bool Porosity_base::readData(std::string userInput) {
        bool tf;
        if (userInput.empty()){
            porInterpolation = interpType::SCALAR;
            porValue = 1.0;
            tf = true;
        }
        else if (ICHNOS::is_input_scalar(userInput)){
            porInterpolation = interpType::SCALAR;
            porValue = std::stod(userInput);
            tf = true;
        }
        else{
            tf = readFile(userInput);
        }
        return tf;
    }

    double Porosity_base::interpolate(vec3 p) {
        double out = 1.0;
        switch (porInterpolation){
            case interpType::SCALAR:
            {
                out = porValue;
                break;
            }
            case interpType::CLOUD:
            {
                out = cloudInterp(p);
                break;
            }
        }
        return out;
    }

    bool Porosity_base::readFile(std::string filename) {
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file" << filename << std::endl;
            return false;
        }
        else{
            std::string line;
            { // Read the TYPE
                getline(datafile, line);
                if (line.compare("CLOUD") == 0){
                    porInterpolation = interpType::CLOUD;
                }
                else{
                    std::cout << "[" << line << "] is unknown interpolation type for porosity" << std::endl;
                    return false;
                }
            }

            {//Read data
                getline(datafile, line);
                std::istringstream inp(line.c_str());
                double rad, pw, thr;
                inp >> rad;
                inp >> pw;
                inp >> thr;
                Radius = rad*rad;
                Power = pw;
                Threshold = thr;

                double x, y, z, v, diam, ratio;
                Pnt_info td;
                int cnt = 0;
                std::vector<cgal_point_3> pp;
                std::vector<Pnt_info> dd;

                double meanRatio = 0.0;
                while (getline(datafile, line)){
                    if (line.size() > 1){
                        std::istringstream inp1(line.c_str());
                        inp1 >> x;
                        inp1 >> y;
                        inp1 >> z;
                        inp1 >> v;
                        inp1 >> td.diameter;
                        inp1 >> td.ratio;
                        td.id = cnt;
                        cnt++;
                        meanRatio += td.ratio;

                        PorosityArray.push_back(v);
                        pp.push_back(cgal_point_3(x, y, z));
                        dd.push_back(td);
                    }
                }
                ratioInit = meanRatio/static_cast<double>(cnt);

                datafile.close();
                auto start = std::chrono::high_resolution_clock::now();
                Tree.insert(boost::make_zip_iterator(boost::make_tuple( pp.begin(),dd.begin() )),
                            boost::make_zip_iterator(boost::make_tuple( pp.end(),dd.end() ) )  );
                Tree.build();
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = finish - start;
                std::cout << "\tPorosity Building time: " << elapsed.count() << std::endl;
                return true;
            }
        }
    }

    double Porosity_base::cloudInterp(vec3 p) {
        cgal_point_3 center(p.x, p.y, p.z);
        std::vector<boost::tuples::tuple<cgal_point_3, Pnt_info>> tmp;
        while(true){
            tmp.clear();
            Fuzzy_sphere fs(center, searchDiameter, searchDiameter/100.0);
            Tree.search(std::back_inserter(tmp), fs);
            if (tmp.size() >= 3){
                break;
            }
            else{
                searchDiameter = searchDiameter*1.5;
                if (searchDiameter > 2*Radius){
                    std::cout << "I have found " << tmp.size() << " points around ("
                              << p.x << "," << p.y << "," << p.z
                              <<") within twice the initial diameter of " << Radius
                              << ". Consider increasing the Initial diameter for Porosity" << std::endl;
                    if (tmp.size() == 0){
                        std::cout << "No Porosity info is found> I will continue with porosity 0.1" << std::endl;
                        return 0.1;
                    }
                    else{
                        break;
                    }
                }
            }
        }

        double xx, yy, zz, actual_dist, scaled_dist, w;
        double sumW = 0;
        double sumWV = 0;
        double mindist = 1000000000;
        double tmp_diam;
        double tmp_ratio;
        for (unsigned int i = 0; i < tmp.size(); ++i){
            xx = p.x - tmp[i].get<0>().x();
            yy = p.y - tmp[i].get<0>().y();
            zz = p.z - tmp[i].get<0>().z();
            actual_dist = std::sqrt(xx * xx + yy * yy + zz * zz);
            if (actual_dist < mindist){
                mindist = actual_dist;
                tmp_diam = tmp[i].get<1>().diameter;
                tmp_ratio = tmp[i].get<1>().ratio;
            }
            zz = zz * ratioCurrent;
            scaled_dist = std::sqrt(xx * xx + yy * yy + zz * zz);
            if (actual_dist < Threshold){
                return PorosityArray[tmp[i].get<1>().id];
            }
            else{
                w = 1 / std::pow(scaled_dist, Power);
                sumW += w;
                sumWV = w*PorosityArray[tmp[i].get<1>().id];
            }
        }
        ratioCurrent = tmp_ratio;
        searchDiameter = tmp_diam;
        return sumWV/sumW;
    }

    void Porosity_base::reset() {
        searchDiameter = Radius;
        ratioCurrent = ratioInit;
    }




}

#endif //ICHNOS_ICHNOS_POROSITY_H
