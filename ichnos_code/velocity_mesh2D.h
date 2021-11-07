//
// Created by giorg on 11/5/2021.
//

#ifndef ICHNOS_VELOCITY_MESH2D_H
#define ICHNOS_VELOCITY_MESH2D_H

#include "ichnos_structures.h"
#include "velocity_base.h"

#if _USEHF > 0
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#endif

namespace ic = ICHNOS;

namespace ICHNOS{
    class Mesh2DVel : public ic::velocityField{
    public:
        Mesh2DVel(boost::mpi::communicator& world_in);

        bool readVelocityField(std::string vf_file, int nPnts);
        void calcVelocity(vec3 &vel,
                          std::vector<int> &ids,
                          std::vector<double> &weights,
                          double time = 0);
        void reset();
        void updateStep(double &step);
        void getVec3Data(std::vector<ic::vec3> &data);
        ic::MeshVelInterpType getInterpType(){return interp_type;};

    private:
        ic::VelTR VEL;
        ic::interpType porType;
        std::vector<std::vector<int>> FaceIds;
        double porosityValue = 1.0;
        ic::MeshVelInterpType interp_type = ic::MeshVelInterpType::UNKNOWN;



        int FrequencyStat;
        int nPoints;
        int nXYZpnts;
        int nSteps;
        int nLayers;
        ic::TimeData tm_data;

        bool readVelocityFiles();
        bool readVH5file(std::string filename, int nPoints, int nSteps);
        bool readFaceIds(std::string filename);
        bool readElementVelocity();
        bool readFaceVelocity();

        void elementInterpolation(vec3& vel,
                                  std::vector<int>& ids,
                                  int i1, int i2, double t);
        void nodeInterpolation(vec3& vel,
                               std::vector<int>& ids,
                               std::vector<double>& weights,
                               double time = 0);
        void faceInterpolation(vec3& vel,
                               std::vector<int>& ids,
                               std::vector<double>& weights,
                               double time = 0);
    };

    Mesh2DVel::Mesh2DVel(boost::mpi::communicator &world_in)
            :
            velocityField(world_in)
    {
        InterpolateOutsideDomain = false;
    }

    bool Mesh2DVel::readVelocityField(std::string vf_file, int nPnts) {
        nXYZpnts = nPnts;
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
                ("Velocity.INTERP", po::value<std::string>(), "Type of interpolation. (ELEMENT, NODE or FACE)")
                ("Velocity.Nlayers", po::value<int>()->default_value(4), "e.g 0002->4, 000->3")
                ("Velocity.TimeStepFile", po::value<std::string>(), "This filename with the time steps")
                ("Velocity.TimeInterp", po::value<std::string>(), "Interpolation type between time steps")
                ("Velocity.RepeatTime", po::value<double>()->default_value(0.0), "The number of days to repeat after the and of time steps")
                ("Velocity.Multiplier", po::value<double>()->default_value(1.0), "This is a multiplier to scale velocity")
                ("Velocity.FaceIdFile", po::value<std::string>(), "Face ids for each element")
                ("MESH2D.Nlayers", po::value<int>()->default_value(4), "Number of layers")
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
            Vtype = ic::castVelType2Enum(vm_vfo["Velocity.Type"].as<std::string>());
            if (Vtype == ic::VelType::INVALID) {
                std::cout << vm_vfo["Velocity.Type"].as<std::string>() << " is an invalid velocity type" << std::endl;
                return false;
            }
            {// Identify the velocity interpolation type
                std::string tp = vm_vfo["Velocity.INTERP"].as<std::string>();
                if (tp.compare("ELEMENT") == 0){
                    interp_type = ic::MeshVelInterpType::ELEMENT;
                }
                else if (tp.compare("NODE") == 0){
                    interp_type = ic::MeshVelInterpType::NODE;
                }
                else if (tp.compare("FACE") == 0){
                    interp_type = ic::MeshVelInterpType::FACE;
                }
                else{
                    std::cout << "The [" << tp << "] is unknown interpolation type" << std::endl;
                    return false;
                }
            }


            multiplier = vm_vfo["Velocity.Multiplier"].as<double>();
            nLayers = vm_vfo["MESH2D.Nlayers"].as<int>();

            // Read the time step
            std::vector<double> TimeSteps;
            if (Vtype == ic::VelType::TRANS){
                std::string TSfile = vm_vfo["Velocity.TimeStepFile"].as<std::string>();
                bool tf = ic::READ::readTimeStepFile(TSfile, TimeSteps);
                if (!tf){return false;}
                std::string TimeInterpType = vm_vfo["Velocity.TimeInterp"].as<std::string>();
                if (TimeInterpType.compare("LINEAR") == 0) {
                    VEL.setTimeInterpolationType(ic::TimeInterpType::LINEAR);
                }
                else{
                    VEL.setTimeInterpolationType(ic::TimeInterpType::NEAREST);
                }
                VEL.setNrepeatDays(vm_vfo["Velocity.RepeatTime"].as<double>());
            }
            else{
                TimeSteps.push_back(0);
                VEL.setTimeInterpolationType(ic::TimeInterpType::NEAREST);
                VEL.setNrepeatDays(0);
            }


            nSteps = TimeSteps.size();

            Prefix = vm_vfo["Velocity.Prefix"].as<std::string>();
            std::string suffix;
            if (vm_vfo.count("Velocity.Suffix")) {
                Suffix = vm_vfo["Velocity.Suffix"].as<std::string>();
            }
            else {
                Suffix = ".ich";
            }
            leadingZeros = vm_vfo["Velocity.LeadingZeros"].as<int>();

            nSteps = TimeSteps.size();

            bool tf = readVelocityFiles();
            if (!tf){return false;}
            VEL.setTSvalue(TimeSteps);

            if (interp_type != ic::MeshVelInterpType::ELEMENT) {
                std::string faceidsfile = vm_vfo["Velocity.FaceIdFile"].as<std::string>();
                tf = readFaceIds(faceidsfile);
                if (!tf) { return false;}
            }
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

    bool Mesh2DVel::readVelocityFiles() {
        switch (interp_type) {
            case ic::MeshVelInterpType::ELEMENT:
            {
                bool tf = readElementVelocity();
                if (!tf){return false;}
                break;
            }
            case ic::MeshVelInterpType::FACE:
            {
                bool tf = readFaceVelocity();
                if (!tf){return false;}
                break;
            }
        }

//#if _USEHF > 0
//        if (suffix.compare(".h5") == 0){
//            std::string fileVXYZ = prefix + ic::num2Padstr(world.rank(), ld_zero) + suffix;
//            bool tf1 = readVH5file(fileVXYZ, nPoints, nSteps);
//            return tf1;
//        }
//#endif
//        std::string fileVFace = prefix + ic::num2Padstr(world.rank(), ld_zero) + suffix;
//        std::vector<std::vector<double>> VFACE;
//        bool tf = ic::READ::read2Darray(fileVFace,nSteps,VFACE);
//        for (int i = 0; i < VFACE[0].size(); i++){
//            for (int j = 0; j < VFACE.size(); ++j){
//                VEL.setVELvalue(VFACE[0][i]*multiplier, i, j, ic::coordDim::vx);
//            }
//        }
        return true;
    }

    bool Mesh2DVel::readVH5file(std::string filename, int nPoints, int nSteps) {
        std::cout << "\tReading file " + filename << std::endl;
#if _USEHF > 0
        const std::string VXYZNameSet("VFACE");
        HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
        HighFive::DataSet datasetVXYZ = HDFNfile.getDataSet(VXYZNameSet);
        std::vector<std::vector<double>> VFACE;
        datasetVXYZ.read(VFACE);
        if (VFACE[0].size() != nPoints){
            std::cout << "The number of points (" << VFACE[0].size() << ") in the file " << filename << " is not " << nPoints << std::endl;
            return false;
        }
        for (int i = 0; i < VFACE[0].size(); i++){
            for (int j = 0; j < VFACE.size(); ++j){
                VEL.setVELvalue(VFACE[j][i]*multiplier, i, j, ic::coordDim::vx);
            }
        }
        return true;
#endif
        return false;
    }

    bool Mesh2DVel::readFaceIds(std::string filename) {
#if _USEHF > 0
        std::string ext = ic::getExtension(filename);
        if (ext.compare(".h5") == 0){
            const std::string FIDNameSet("FACEID");
            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet datasetFID = HDFNfile.getDataSet(FIDNameSet);
            std::vector<std::vector<int>> fids;
            datasetFID.read(fids);
            for (int i = 0; i < fids[0].size(); i++){
                std::vector<int> tmp_ids;
                for (int j = 0; j < fids.size(); ++j){
                    if (fids[j][i] != 0)
                        tmp_ids.push_back(fids[j][i]);
                }
                FaceIds.push_back(tmp_ids);
            }
            return true;
        }
#endif
        {
            std::vector<std::vector<int>> fids;
            bool tf = ic::READ::read2Darray(filename,4,fids);
            if (!tf){return false;}
            for (int i = 0; i < fids.size(); i++){
                std::vector<int> tmp_ids;
                for (int j = 0; j < fids[i].size(); ++j){
                    if (fids[i][j] != 0)
                        tmp_ids.push_back(fids[i][j]);
                }
                FaceIds.push_back(tmp_ids);
            }
            return true;
        }
    }

    void Mesh2DVel::calcVelocity(vec3 &vel, std::vector<int> &ids, std::vector<double> &weights, double tm) {
        if (ids[0] == -9 || ids[1] == -9){
            vel = 99999;
            return;
        }
        int i1, i2;
        double t;
        VEL.findIIT(tm, i1, i2, t);
        switch (interp_type){
            case ic::MeshVelInterpType::ELEMENT:
            {
                elementInterpolation(vel, ids, i1, i2, t);
                break;
            }
            case ic::MeshVelInterpType::NODE:
            {
                nodeInterpolation(vel, ids, weights, tm);
                break;
            }
            case ic::MeshVelInterpType::FACE:
            {
                faceInterpolation(vel, ids, weights, tm);
                break;
            }
        }

        double porosity = 1.0;
        if (porType == ic::interpType::CLOUD) {
            // TODO porosity = ic::interpolateScalarTree(PorosityTree, p);
        }
        else if (porType == ic::interpType::SCALAR)
            porosity = porosityValue;
        vel = vel * (1/porosity);

        tm_data.tm = VEL.getTSvalue(i1) * (1-t) + VEL.getTSvalue(i2)*t;
        tm_data.idx1 = i1;
        tm_data.idx2 = i2;
        tm_data.t = t;
    }

    void Mesh2DVel::updateStep(double& step) {

    }

    void Mesh2DVel::getVec3Data(std::vector<ic::vec3>& data) {

    }

    void Mesh2DVel::reset() {

    }

    void Mesh2DVel::elementInterpolation(vec3 &vel,
                                         std::vector<int> &ids,
                                         int i1, int i2, double t) {
        // In element interpolation we only need the ids
        int elid = ids[0];
        int lay = ids[1];
        int idx = elid + lay*nXYZpnts;
        vel = VEL.getVelocity(idx, i1, i2, t);
    }

    void Mesh2DVel::nodeInterpolation(vec3& vel,
                                      std::vector<int>& ids,
                                      std::vector<double>& weights,
                                      double time) {
        int elid = ids[0];
        int lay = ids[1];


    }
    void Mesh2DVel::faceInterpolation(ic::vec3& vel,
                                      std::vector<int>& ids,
                                      std::vector<double>& weights,
                                      double time) {
        int elid = ids[0];
        int lay = ids[1];


    }

    bool Mesh2DVel::readElementVelocity() {
        switch (Vtype) {
            case ic::VelType::STEADY:
            {
                std::string filename = Prefix + ic::num2Padstr(world.rank(), leadingZeros) + Suffix;
                if (Suffix.compare(".h5") == 0){
                    std::cout << "MESH2D - Steady State - H5 input is not implemented yet" << std::endl;
                    return false;
                }
                std::cout << "\tReading file " + filename << std::endl;
                std::ifstream datafile(filename.c_str());
                if (!datafile.good()){
                    //TODO
                    std::cout << "Can't open the file " << filename << std::endl;
                    return false;
                }
                else{
                    std::string line;
                    double vx, vy, vz;
                    std::vector<ic::vec3> tmp_vel;
                    while (getline(datafile, line)){
                        if (line.size() > 1){
                            std::istringstream inp(line.c_str());
                            inp >> vx;
                            inp >> vy;
                            inp >> vz;
                            tmp_vel.push_back(ic::vec3(vx, vy, vz));
                        }
                    }
                    datafile.close();
                    nPoints = tmp_vel.size();
                    VEL.init(nPoints, nSteps);
                    for (unsigned int i = 0; i < tmp_vel.size(); ++i){
                        VEL.setVELvalue(tmp_vel[i].x * multiplier, i, 0, ic::coordDim::vx);
                        VEL.setVELvalue(tmp_vel[i].y * multiplier, i, 0, ic::coordDim::vy);
                        VEL.setVELvalue(tmp_vel[i].z * multiplier, i, 0, ic::coordDim::vz);
                    }
                    return true;
                }
                break;
            }
            case ic::VelType::TRANS:
            {
                //TODO
                break;
            }

        }



        return true;
    }

    bool Mesh2DVel::readFaceVelocity() {
        return true;
    }

}

#endif //ICHNOS_VELOCITY_MESH2D_H
