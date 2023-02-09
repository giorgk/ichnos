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
        double value(vec3 p);

    private:
        interpType porInterpolation;
        double porValue;
        bool readFile(std::string filename);
        search_tree_info Tree;
        std::vector<double> PorosityArray;

        double Radius;
        double Power;
    };

    bool Porosity_base::readData(std::string userInput) {
        if (ICHNOS::is_input_scalar(userInput)){
            porInterpolation = interpType::SCALAR;
            porValue = std::stod(userInput);
        }
        else{
            bool tf = readFile(userInput);

        }
    }

    double Porosity_base::value(vec3 p) {
        double out = 1.0;
        switch (porInterpolation){
            case interpType::SCALAR:
            {
                out = porValue;
                break;
            }
            case interpType::CLOUD:
            {


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
                double rad, pw;
                inp >> rad;
                inp >> pw;
                Radius = rad*rad;
                Power = pw;

                double x, y, z, v, diam, ratio;
                Pnt_info td;
                int cnt = 0;
                std::vector<cgal_point_3> pp;
                std::vector<Pnt_info> dd;

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

                        PorosityArray.push_back(v);
                        pp.push_back(cgal_point_3(x, y, z));
                        dd.push_back(td);
                    }
                }
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




}

#endif //ICHNOS_ICHNOS_POROSITY_H
