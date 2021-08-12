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

    struct TimeData{
        double tm;
        int idx1;
        int idx2;
        double t;
        double dir;
    };

    class transVel : public ic::velocityField{
    public:
        transVel(boost::mpi::communicator& world_in);
        void readVelocityField(std::string vf_file, int nPnts);
        void calcVelocity(ic::vec3& vel,
                          std::vector<int>& ids,
                          std::vector<double>& weights,
                          double time = 0);
        void reset();
        void updateStep(double& step);
        void getVec3Data(std::vector<ic::vec3>& data);

    private:
        //bool readXYZfile(std::string prefix, std::string suffix, int ld_zero);
        bool readVelocityFiles(std::string prefix, std::string suffix, int ld_zero);
        bool readVfile(std::string filename, int nPoints, int nSteps, ic::coordDim dim);
        bool readSteadyVfile(std::string filename, int nPoints);
        bool readTimeSteps(std::string filename, std::vector<double>& TS);

        ic::VelTR VEL;
        ic::search_tree_info Tree;

        ic::interpType porType;
        double porosityValue = 1.0;


        int nPoints;
        int nSteps;

        double multiplier = 1.0;
        //double Power;
        //double Scale = 1.0;
        //double initial_diameter = 640;
        //double initial_ratio = 20;
        //double diameter;
        //double ratio;

        //double Threshold;
        int FrequencyStat;
        double calc_time = 0.0;
        double max_calc_time = 0.0;
        int count_times = 0;

        TimeData tm_data;

        bool bIsInitialized = false;
        //double search_mult = 2.5;
        ic::vec3 ll, uu, pp, vv;
        ic::TimeInterpType timeInterp = ic::TimeInterpType::NEAREST;

    };

    transVel::transVel(boost::mpi::communicator &world_in)
        :
        velocityField(world_in)
    {
        InterpolateOutsideDomain = true;
    }

    void transVel::readVelocityField(std::string vf_file, int nPnts) {
        if (world.rank() == 0)
            std::cout << "--> Velocity configuration file: " << vf_file << std::endl;

        po::options_description velocityFieldOptions("Velocity field options");
        po::variables_map vm_vfo;
        velocityFieldOptions.add_options()
            // Velocity parameters
            ("Velocity.Prefix", po::value<std::string>(), "Prefix for the filename")
            ("Velocity.LeadingZeros", po::value<int>()->default_value(4), "e.g 0002->4, 000->3")
            ("Velocity.Suffix", po::value<std::string>(), "ending of file after procid")
            ("Velocity.Type", po::value<std::string>(), "Type of velocity. (STEADY or TRANS)")
            ("Velocity.TimeStepFile", po::value<std::string>(), "This filename with the time steps")
            ("Velocity.TimeInterp", po::value<std::string>(), "Interpolation type between timesteps")
            ("Velocity.Multiplier", po::value<double>()->default_value(1), "This is a multiplier to scale velocity")
            //("Velocity.Scale", po::value<double>()->default_value(1.0), "Scale the domain before velocity calculation")
            //("Velocity.Power", po::value<double>()->default_value(3.0), "Power of the IDW interpolation")
            //("Velocity.InitDiameter", po::value<double>()->default_value(5000), "Initial diameter")
            //("Velocity.InitRatio", po::value<double>()->default_value(1), "Initial ratio")



            // Porosity parameters
            ("Porosity.Value", po::value<std::string>(), "Porosity. Either a file or a single number")

            //General
            ("General.OwnerThreshold", po::value<double>()->default_value(0.75), "Threshold for the processor ownership")
            //("General.Threshold", po::value<double>()->default_value(0.001), "Threshold of distance of IDW")
            ("General.FrequencyStat", po::value<int>()->default_value(20), "Frequency of printing stats")
        ;

        po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions, true), vm_vfo);

        { // General
            OwnerThreshold = vm_vfo["General.OwnerThreshold"].as<double>();
            //Threshold = vm_vfo["General.Threshold"].as<double>();
            FrequencyStat = vm_vfo["General.FrequencyStat"].as<int>();
        }

        {// Velocity parameters
            multiplier = vm_vfo["Velocity.Multiplier"].as<double>();
            //Scale = vm_vfo["Velocity.Scale"].as<double>();
            //Power = vm_vfo["Velocity.Power"].as<double>();
            //initial_diameter = vm_vfo["Velocity.InitDiameter"].as<double>();
            //initial_ratio = vm_vfo["Velocity.InitRatio"].as<double>();

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

            Vtype = ic::castVelType2Enum(vm_vfo["Velocity.Type"].as<std::string>());
            if (Vtype == ic::VelType::INVALID) {
                std::cout << vm_vfo["Velocity.Type"].as<std::string>() << " is an invalid velocity type" << std::endl;
                return;
            }


            nPoints = nPnts;
            VEL.init(nPoints, nSteps);
            VEL.setTSvalue(TimeSteps);

            // Read the Velocity files
            tf = readVelocityFiles(prefix, suffix, leadZeros);
        }

        { // Porosity
            if (vm_vfo.count("Porosity.Value")){
                std::string porfile = vm_vfo["Porosity.Value"].as<std::string>();
                if (porfile.empty()){
                    porType = ic::interpType::INGORE;
                }
                else{
                    if (ic::is_input_scalar(porfile)) {
                        porType = ic::interpType::SCALAR;
                        porosityValue = std::stod(porfile);
                    }
                    else{
                        // TODO
                    }
                }
            }

        }
    }

    bool transVel::readVelocityFiles(std::string prefix, std::string suffix, int ld_zero){
        bool tf = false;
        auto start = std::chrono::high_resolution_clock::now();
        if (Vtype == ic::VelType::TRANS){
            std::string fileVX = prefix + "VX_" + ic::num2Padstr(/*dbg_rank*/world.rank(), ld_zero) + suffix;
            tf = readVfile(fileVX, nPoints, nSteps, ic::coordDim::vx);
            if (tf){
                std::string fileVY = prefix + "VY_" + ic::num2Padstr(/*dbg_rank*/world.rank(), ld_zero) + suffix;
                tf = readVfile(fileVY, nPoints, nSteps, ic::coordDim::vy);
            }
            if (tf){
                std::string fileVZ = prefix + "VZ_" + ic::num2Padstr(/*dbg_rank*/world.rank(), ld_zero) + suffix;
                tf = readVfile(fileVZ, nPoints, nSteps, ic::coordDim::vz);
            }
        }
        else if (Vtype == ic::VelType::STEADY){
            std::string fileVX = prefix + ic::num2Padstr(world.rank(), ld_zero) + suffix;
            tf = readSteadyVfile(fileVX,nPoints);
        }

        if (!tf)
            return false;

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "\tRead Velocity in : " << elapsed.count() << std::endl;
        return true;
    }

    bool transVel::readVfile(std::string filename, int nPoints, int nSteps, ic::coordDim dim){
        std::cout << "\tReading file " + filename << std::endl;
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

    bool transVel::readSteadyVfile(std::string filename, int nPoints){
        std::cout << "Reading file " + filename << std::endl;
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()){
            std::cout << "Can't open the file " << filename << std::endl;
            return false;
        }
        else{
            std::string line;
            double coord, vx, vy, vz;
            int tmp;
            for (int i = 0; i < nPoints; i++){
                getline(datafile, line);
                std::istringstream inp(line.c_str());
                inp >> coord;
                inp >> coord;
                inp >> coord;
                inp >> tmp;
                inp >> coord;
                inp >> coord;
                inp >> vx;
                inp >> vy;
                inp >> vz;
                VEL.setVELvalue(vx*multiplier, i, 0, ic::coordDim::vx);
                VEL.setVELvalue(vy*multiplier, i, 0, ic::coordDim::vy);
                VEL.setVELvalue(vz*multiplier, i, 0, ic::coordDim::vz);

            }
        }
        return true;
    }

    bool transVel::readTimeSteps(std::string filename, std::vector<double>& TS){
        std::cout << "\tReading file " + filename << std::endl;
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
            std::cout << "\tProc["<< world.rank() <<  "] Number of time steps: " << nSteps << std::endl;
            datafile.close();
        }
        return true;
    }

//    bool transVel::readXYZfile(std::string prefix, std::string suffix, int ld_zero){
//        std::string fileXYZ = prefix + "XYZ_" + ic::num2Padstr(/*dbg_rank*/world.rank(), ld_zero) + suffix;
//        std::cout << "Reading file " + fileXYZ << std::endl;
//        std::ifstream datafile(fileXYZ.c_str());
//        if (!datafile.good()){
//            std::cout << "Can't open the file " << fileXYZ << std::endl;
//            return false;
//        }
//        else{
//            std::vector<ic::cgal_point_3> pp;
//            std::vector<ic::Pnt_info> dd;
//
//            std::string line;
//            double x, y, z;
//            ic::Pnt_info td;
//            int cnt = 0;
//            while (getline(datafile, line)){
//                if (line.size() > 1){
//                    std::istringstream inp(line.c_str());
//                    inp >> x;
//                    inp >> y;
//                    inp >> z;
//                    inp >> td.proc;
//                    inp >> td.diameter;
//                    inp >> td.ratio;
//                    td.id = cnt;
//                    cnt++;
//
//                    pp.push_back(ic::cgal_point_3(x, y, z));
//                    dd.push_back(td);
//                }
//            }
//            nPoints = pp.size();
//            std::cout << "Proc["<< world.rank() <<  "] Number of XYZ points: " << nPoints << std::endl;
//            datafile.close();
//
//            { //Build tree
//                auto start = std::chrono::high_resolution_clock::now();
//                Tree.insert(boost::make_zip_iterator(boost::make_tuple( pp.begin(),dd.begin() )),
//                            boost::make_zip_iterator(boost::make_tuple( pp.end(),dd.end() ) )  );
//                Tree.build();
//                auto finish = std::chrono::high_resolution_clock::now();
//                std::chrono::duration<double> elapsed = finish - start;
//                std::cout << "Point Set Building time: " << elapsed.count() << std::endl;
//            }
//        }
//        return true;
//    }

    void transVel::calcVelocity(ic::vec3& vel,
                                std::vector<int>& ids,
                                std::vector<double>& weights,
                                double tm){


        // Find the time step
        int i1, i2;
        double t;
        VEL.findIIT(tm, i1, i2, t);

        double sumW = 0;
        ic::vec3 sumWVal;
        vel.zero();
        for (unsigned int i = 0; i < ids.size(); ++i) {
            vel = VEL.getVelocity(ids[i], i1, i2, t, timeInterp);
            sumW += weights[i];
            sumWVal = sumWVal + vel * weights[i];
        }
        vel = sumWVal * (1 / sumW);

        {
//            auto start = std::chrono::high_resolution_clock::now();
//
//
//            auto finish = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double> elapsed = finish - start;
//            calc_time += elapsed.count();
//            if (elapsed.count() > max_calc_time)
//                max_calc_time = elapsed.count();
        }

        double porosity = 1.0;
        if (porType == ic::interpType::CLOUD) {
            // TODO porosity = ic::interpolateScalarTree(PorosityTree, p);
        }
        else if (porType == ic::interpType::SCALAR)
            porosity = porosityValue;
        vel = vel * (1/porosity);

        // TODO Understand why this is needed
        vv = vel;
        tm_data.tm = tm;
        tm_data.idx1 = i1;
        tm_data.idx2 = i2;
        tm_data.t = t;

//        count_times++;
//        ic::PrintStat(count_times, FrequencyStat, calc_time, max_calc_time);
    }

    void transVel::reset() {
        bIsInitialized = false;
    }

    void transVel::getVec3Data(std::vector<ic::vec3> &data) {
        pp = data[0];
        ll = data[1];
        uu = data[2];
    }

    void transVel::updateStep(double &step) {
        double stepLen, stepTime;
        double dst = ic::diameter_along_velocity(pp, vv, ll, uu);

        if (stepOpt.nSteps > 0){
            // This is the length that respects the nSteps
            stepLen = dst/stepOpt.nSteps;
        }
        else{
            stepLen = 10000000;
        }

        if (stepOpt.StepSize > 0){
            // If the length of the user step is smaller than the step defined by the nSteps.
            // Use the used defined that is smaller
            if (stepOpt.StepSize < stepLen){
                stepLen = stepOpt.StepSize;
            }
        }

        if (stepOpt.StepSizeTime > 0){
            stepTime = vv.len()*stepOpt.StepSizeTime;
        }
        else{
            stepTime = 10000000;
        }

        if (stepOpt.nStepsTime > 0){
            if (tm_data.idx1 != tm_data.idx2){
                double dt = 1.0/ static_cast<double>(stepOpt.nStepsTime);
                double tm_1 = VEL.getTSvalue(tm_data.idx1);
                double tm_2 = VEL.getTSvalue(tm_data.idx2);
                double tmp_step = dt*(tm_2 - tm_1);
                double end_time = tm_data.tm + stepOpt.dir*tmp_step;
                if (stepOpt.dir > 0){
                    if (end_time > tm_2 && tm_data.idx2 < nSteps - 1){
                        stepTime = vv.len() * (tm_2 - tm_data.tm);
                    }
                    else{
                        stepTime = vv.len() * tmp_step;
                    }
                }
                else{
                    if (end_time < tm_1 && tm_data.idx1 > 0){
                        stepTime = vv.len() * (tm_data.tm - tm_1);
                    }
                    else{
                        stepTime = vv.len() * tmp_step;
                    }
                }
            }
        }
        step = std::min<double>(stepTime, stepLen);
    }
}

