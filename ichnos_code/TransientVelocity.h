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

#if _USEHF > 0
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#endif

namespace po = boost::program_options;
namespace ic = ICHNOS;

namespace ICHNOS{
    class CloudVel : public ic::velocityField{
    public:
        CloudVel(boost::mpi::communicator& world_in, XYZ_base &XYZ_in);
        bool readVelocityField(std::string vf_file);
        void calcVelocity(vec3& p, vec3& vel,
                          std::map<int, double>& proc_map,
                          helpVars& pvlu,
                          bool& out,
                          //std::vector<int>& ids,
                          //std::vector<double>& weights,
                          double time = 0);
        void reset(Streamline& S);
        double stepTimeupdate(helpVars& pvlu);
        void updateStep(double& step);
        void getVec3Data(std::vector<ic::vec3>& data);

    protected:
        //bool readXYZfile(std::string prefix, std::string suffix, int ld_zero);
        bool readVelocityFiles();
        void setVelData(std::vector<std::vector<double>>& data, double multiplier, coordDim dim);
        //bool readVfile(std::string filename, int nPoints, int nSteps, ic::coordDim dim);
        //bool readVH5file(std::string filename, int nPoints, int nSteps);
        //bool readSteadyVfile(std::string filename, int nPoints);
        //bool readTimeSteps(std::string filename, std::vector<double>& TS);


        ic::VelTR VEL;
        //ic::search_tree_info Tree;

        ic::interpType porType;
        double porosityValue = 1.0;


        int nPoints;
        int nSteps;

        //double Power;
        //double Scale = 1.0;
        //double initial_diameter = 640;
        //double initial_ratio = 20;
        //double diameter;
        //double ratio;

        //double Threshold;
        //int FrequencyStat;
        //double calc_time = 0.0;
        //double max_calc_time = 0.0;
        //int count_times = 0;

        //ic::TimeData tm_data;

        //bool bIsInitialized = false;
        //double search_mult = 2.5;
        //ic::vec3 ll, uu, pp, vv;
        ic::TimeInterpType timeInterp = ic::TimeInterpType::NEAREST;

    };

    CloudVel::CloudVel(boost::mpi::communicator &world_in, XYZ_base &XYZ_in)
        :
        velocityField(world_in, XYZ_in)
    {
        InterpolateOutsideDomain = true;
    }

    bool CloudVel::readVelocityField(std::string vf_file) {
        if (world.rank() == 0)
            std::cout << "--> Velocity configuration file: " << vf_file << std::endl;

        po::options_description velocityFieldOptions("Velocity field options");
        po::variables_map vm_vfo;
        velocityFieldOptions.add_options()
            // Velocity parameters
            ("Velocity.Prefix", po::value<std::string>(), "Prefix for the filename")
            ("Velocity.LeadingZeros", po::value<int>()->default_value(4), "e.g 0002->4, 000->3")
            ("Velocity.Suffix", po::value<std::string>(), "ending of file after procid")
            ("Velocity.Type", po::value<std::string>(), "Type of velocity.")
            //("Velocity.Trans", po::value<int>()->default_value(0), "0->steady state, 1->Transient state")
            ("Velocity.TimeStepFile", po::value<std::string>(), "This filename with the time steps")
            ("Velocity.TimeInterp", po::value<std::string>(), "Interpolation type between time steps")
            ("Velocity.RepeatTime", po::value<double>()->default_value(0.0), "The number of days to repeat after the end of time steps")
            ("Velocity.Multiplier", po::value<double>()->default_value(1.0), "This is a multiplier to scale velocity")

            // Porosity parameters
            ("Porosity.Value", po::value<std::string>(), "Porosity. Either a file or a single number")

            //General
            ("Other.OwnerThreshold", po::value<double>()->default_value(0.75), "Threshold for the processor ownership")
            //("General.Threshold", po::value<double>()->default_value(0.001), "Threshold of distance of IDW")
            ("Other.FrequencyStat", po::value<int>()->default_value(20), "Frequency of printing stats")
        ;

        po::store(po::parse_config_file<char>(vf_file.c_str(), velocityFieldOptions, true), vm_vfo);

        { // General
            OwnerThreshold = vm_vfo["Other.OwnerThreshold"].as<double>();
            //Threshold = vm_vfo["General.Threshold"].as<double>();
            //FrequencyStat = vm_vfo["General.FrequencyStat"].as<int>();
        }

        {// Velocity parameters
            Vtype = ic::castVelType2Enum(vm_vfo["Velocity.Type"].as<std::string>());
            if (Vtype == ic::VelType::INVALID) {
                std::cout << vm_vfo["Velocity.Type"].as<std::string>() << " is an invalid velocity type" << std::endl;
                return false;
            }

            multiplier = vm_vfo["Velocity.Multiplier"].as<double>();
            //nPoints = vm_vfo["Velocity.nPoints"].as<int>();


            //Scale = vm_vfo["Velocity.Scale"].as<double>();
            //Power = vm_vfo["Velocity.Power"].as<double>();
            //initial_diameter = vm_vfo["Velocity.InitDiameter"].as<double>();
            //initial_ratio = vm_vfo["Velocity.InitRatio"].as<double>();

            //isVeltrans  = vm_vfo["Velocity.Trans"].as<int>() != 0;
            // Read the time step
            std::vector<double> TimeSteps;

            isVeltrans = false;//  = vm_vfo["Velocity.Trans"].as<int>() != 0;
            if (vm_vfo.count("Velocity.TimeStepFile")){
                std::string TSfile = vm_vfo["Velocity.TimeStepFile"].as<std::string>();
                if (!TSfile.empty()){
                    isVeltrans = true;
                }
            }

            if (isVeltrans){
                std::string TSfile = vm_vfo["Velocity.TimeStepFile"].as<std::string>();
                bool tf = ic::READ::readTimeStepFile(TSfile, TimeSteps);
                if (!tf) {return false;}
                std::string TimeInterpType = vm_vfo["Velocity.TimeInterp"].as<std::string>();
                if (TimeInterpType.compare("LINEAR") == 0) {
                    timeInterp = ic::TimeInterpType::LINEAR;
                }
                VEL.setTimeInterpolationType(timeInterp);
                VEL.setNrepeatDays(vm_vfo["Velocity.RepeatTime"].as<double>());
            }
            else{
                TimeSteps.push_back(0);
                timeInterp = ic::TimeInterpType::NEAREST;
                VEL.setNrepeatDays(0);
                VEL.setTimeInterpolationType(timeInterp);
            }

            Prefix = vm_vfo["Velocity.Prefix"].as<std::string>();
            if (vm_vfo.count("Velocity.Suffix")) {
                Suffix = vm_vfo["Velocity.Suffix"].as<std::string>();
            }
            else {
                Suffix = ".ich";
            }
            leadingZeros = vm_vfo["Velocity.LeadingZeros"].as<int>();


            nSteps = TimeSteps.size()-1;
            VEL.setnSteps(nSteps);
            //if (nPoints > 0){
            //    // We can initialize at this point only if the nPoints is given.
            //    // Otherwise, we will initialize it inside the readers
            //    VEL.init(nPoints);
            //}
            VEL.setTSvalue(TimeSteps);

            // Read the Velocity files
            bool tf = readVelocityFiles();
            if (!tf)
                return false;
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
        return true;
    }

    bool CloudVel::readVelocityFiles(){
        bool tf = false;
        int proc_id = world.rank();
        if (XYZ.runAsThread){
            proc_id = 0;
        }

        auto start = std::chrono::high_resolution_clock::now();

        if (isVeltrans){

            if (Suffix.compare(".h5") == 0){
#if _USEHF > 0
                std::string fileVXYZ = Prefix + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
                tf = READ::H5Transient3DVelocity(fileVXYZ, nSteps, multiplier, VEL);
                nPoints = VEL.getNpoints();
#endif
            }
            else{
                { //VX
                    std::string fileVX = Prefix + "VX_" + ic::num2Padstr(/*dbg_rank*/proc_id, leadingZeros) + Suffix;
                    std::vector<std::vector<double>> data;
                    tf = READ::read2Darray<double>(fileVX, nSteps, data);
                    if (tf) {
                        nPoints = static_cast<int>(data.size());
                        VEL.init(nPoints, nSteps, 3);
                        setVelData(data, multiplier, ic::coordDim::vx);
                    }
                    else{
                        tf = false;
                    }
                }


                if (tf){//VY
                    std::string fileVY = Prefix + "VY_" + ic::num2Padstr(/*dbg_rank*/proc_id, leadingZeros) + Suffix;
                    std::vector<std::vector<double>> data;
                    tf = READ::read2Darray<double>(fileVY, nSteps, data);
                    if (tf){
                        if (nPoints != data.size()){
                            std::cout << "The number of data (" << nPoints << ") in VX file is different than the VY file (" << data.size() << ")" << std::endl;
                            tf = false;
                        }
                        setVelData(data, multiplier, ic::coordDim::vy);
                    }
                    else{
                        tf = false;
                    }
                }

                if (tf){// VZ
                    std::string fileVZ = Prefix + "VZ_" + ic::num2Padstr(/*dbg_rank*/proc_id, leadingZeros) + Suffix;
                    std::vector<std::vector<double>> data;
                    tf = READ::read2Darray<double>(fileVZ, nSteps, data);
                    if (tf){
                        if (nPoints != data.size()){
                            std::cout << "The number of data (" << nPoints << ") in VX and VY file is different than the VZ file (" << data.size() << ")" << std::endl;
                            tf = false;
                        }
                        setVelData(data, multiplier, ic::coordDim::vz);
                        tf = true;
                    }
                    else{
                        tf = false;
                    }
                }
            }
        }
        else {
            std::string fileVXYZ = Prefix + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
            if (Suffix.compare(".h5") == 0){
#if _USEHF > 0
                std::vector<std::vector<double>> VXYZ;
                tf = READ::H5SteadyState3DVelocity(fileVXYZ, VXYZ);
                if (tf){
                    nPoints = VXYZ[0].size();
                    VEL.init(nPoints,nSteps, 3);
                    for (int i = 0; i < VXYZ[0].size(); i++){
                        VEL.setVELvalue(VXYZ[0][i]*multiplier, i, 0, coordDim::vx);
                        VEL.setVELvalue(VXYZ[1][i]*multiplier, i, 0, coordDim::vy);
                        VEL.setVELvalue(VXYZ[2][i]*multiplier, i, 0, coordDim::vz);
                    }
                    tf = true;
                }
#endif
            }
            else{
                std::vector<std::vector<double>> data;
                tf = READ::read2Darray<double>(fileVXYZ, 9, data);
                if (tf){
                    nPoints = static_cast<int>(data.size());
                    VEL.init(nPoints, 1, 3);
                    for (int i = 0; i < nPoints; i++){
                        VEL.setVELvalue(data[i][6]*multiplier, i, 0, coordDim::vx);
                        VEL.setVELvalue(data[i][7]*multiplier, i, 0, coordDim::vy);
                        VEL.setVELvalue(data[i][8]*multiplier, i, 0, coordDim::vz);
                    }
                }
            }
        }

        if (tf){
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            std::cout << "\tRead Velocity in : " << elapsed.count() << std::endl;
        }
        return tf;
    }

    void CloudVel::setVelData(std::vector<std::vector<double>>& data, double multiplier, coordDim dim){
        for (int i = 0; i < nPoints; i++){
            for (int j = 0; j < nSteps; j++){
                VEL.setVELvalue(data[i][j]*multiplier, i, j, dim);
            }
        }
    }

    /*
    bool CloudVel::readVH5file(std::string filename, int nPoints, int nSteps){
        std::cout << "\tReading file " + filename << std::endl;
#if _USEHF > 0
        if (!isVeltrans){
            const std::string VXYZNameSet("VXYZ");
            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet datasetVXYZ = HDFNfile.getDataSet(VXYZNameSet);
            std::vector<std::vector<double>> VXYZ;
            datasetVXYZ.read(VXYZ);
            if (VXYZ[0].size() != nPoints){
                std::cout << "The number of points (" << VXYZ[0].size() << ") in the file " << filename << " is not " << nPoints << std::endl;
                return false;
            }
            for (int i = 0; i < VXYZ[0].size(); i++){
                VEL.setVELvalue(VXYZ[0][i]*multiplier, i, 0, ic::coordDim::vx);
                VEL.setVELvalue(VXYZ[1][i]*multiplier, i, 0, ic::coordDim::vy);
                VEL.setVELvalue(VXYZ[2][i]*multiplier, i, 0, ic::coordDim::vz);
            }
        }
        else {
            const std::string VXNameSet("VX");
            const std::string VYNameSet("VY");
            const std::string VZNameSet("VZ");
            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet datasetVX = HDFNfile.getDataSet(VXNameSet);
            HighFive::DataSet datasetVY = HDFNfile.getDataSet(VYNameSet);
            HighFive::DataSet datasetVZ = HDFNfile.getDataSet(VZNameSet);
            std::vector<std::vector<double>> VX;
            std::vector<std::vector<double>> VY;
            std::vector<std::vector<double>> VZ;
            datasetVX.read(VX);
            datasetVY.read(VY);
            datasetVZ.read(VZ);

            if (VX[0].size() != nPoints){
                std::cout << "The number of points (" << VX[0].size() << ") in the file " << filename << " is not " << nPoints << std::endl;
                return false;
            }
            if (VX.size() != nSteps){
                std::cout << "The number of time steps (" << VX.size() << ") in the file " << filename << " is not " << nSteps << std::endl;
                return false;
            }

            for (int i = 0; i < VX[0].size(); i++){
                for (int j = 0; j < VX.size(); ++j) {
                    VEL.setVELvalue(VX[j][i]*multiplier, i, j, ic::coordDim::vx);
                    VEL.setVELvalue(VY[j][i]*multiplier, i, j, ic::coordDim::vy);
                    VEL.setVELvalue(VZ[j][i]*multiplier, i, j, ic::coordDim::vz);
                }
            }
        }
        return true;
#endif
        return false;
    }
*/
    /*
    bool CloudVel::readVfile(std::string filename, int nPoints, int nSteps, ic::coordDim dim){
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
     */

    /*
    bool CloudVel::readSteadyVfile(std::string filename, int nPoints){
        std::cout << "\tReading file " + filename << std::endl;
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
     */

    /*
    bool CloudVel::readTimeSteps(std::string filename, std::vector<double>& TS){
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
*/

//    bool CloudVel::readXYZfile(std::string prefix, std::string suffix, int ld_zero){
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

    void CloudVel::calcVelocity(vec3& p, vec3& vel,
                                std::map<int, double>& proc_map,
                                helpVars& pvlu,
                                bool& out,
                                //std::vector<int>& ids,
                                //std::vector<double>& weights,
                                double tm){

        std::vector<int> ids;
        std::vector<double> weights;
        XYZ.calcWeights(p, ids,weights, proc_map, pvlu, out);
        if (!out)
            return;


        // Find the time step
        int i1, i2;
        double t, tm_tmp;
        if (isVeltrans){
            VEL.findIIT(tm, i1, i2, t, tm_tmp);
        }
        else{
            i1 = 0;
            i2 = 0;
            t = 0.0;
            tm_tmp = 0.0;
        }
        pvlu.td.idx1 = i1;
        pvlu.td.idx2 = i2;
        pvlu.td.t = t;
        pvlu.td.tm = tm;
        pvlu.td.tm_tmp = tm_tmp;

        double sumW = 0;
        ic::vec3 sumWVal;
        vel.zero();
        for (unsigned int i = 0; i < ids.size(); ++i) {
            vel = VEL.getVelocity(ids[i], i1, i2, t);
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

        //vv = vel;
        //tm_data.tm = VEL.getTSvalue(i1) * (1-t) + VEL.getTSvalue(i2)*t;
        //tm_data.idx1 = i1;
        //tm_data.idx2 = i2;
        //tm_data.t = t;

//        count_times++;
//        ic::PrintStat(count_times, FrequencyStat, calc_time, max_calc_time);
    }

    void CloudVel::reset(Streamline& S) {
        //bIsInitialized = false;
    }

    void CloudVel::getVec3Data(std::vector<ic::vec3> &data) {
        //pp = data[0];
        //ll = data[1];
        //uu = data[2];
    }

    double CloudVel::stepTimeupdate(helpVars &pvlu) {
        return VEL.stepTimeUpdate(pvlu, stepOpt);
    }

    void CloudVel::updateStep(double &step) {
        std::cout << "The Code should never be here!!!" << std::endl;
        /*
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
                    if (std::abs(tm_data.tm - tm_2) < 0.25*tmp_step){
                        if (tm_data.idx2 + 1 < nSteps){
                            tm_2 = VEL.getTSvalue(tm_data.idx2+1);
                        }
                    }
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
        */
    }
}

