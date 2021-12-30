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

        std::vector<vec3> nds;
        std::vector<std::vector<int>> msh;



        int FrequencyStat;
        int nPoints;
        int nSteps;
        int nLayers;
        int nFaces;
        int nElements;
        int nNodes;
        int nTotalFaces;
        ic::TimeData tm_data;

        bool readVelocityFiles();
        bool readVH5file(std::string filename, int nPoints, int nSteps);
        bool readFaceIds(std::string filename);
        bool readXYZVelocity();
        bool readFaceVelocity();

        void elementInterpolation(vec3& vel,
                                  std::vector<int>& ids,
                                  int i1, int i2, double t);
        void nodeInterpolation(vec3& vel,
                               std::vector<int>& ids,
                               std::vector<double>& weights,
                               int i1, int i2, double t);
        void faceInterpolation(vec3& vel,
                               std::vector<int>& ids,
                               std::vector<double>& weights,
                               int i1, int i2, double t);
    };

    Mesh2DVel::Mesh2DVel(boost::mpi::communicator &world_in)
            :
            velocityField(world_in)
    {
        InterpolateOutsideDomain = false;
    }

    bool Mesh2DVel::readVelocityField(std::string vf_file, int nPnts) {

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
                ("MESH2D.FaceIdFile", po::value<std::string>(), "Face ids for each element, Required for FACE type")
                ("MESH2D.Nlayers", po::value<int>()->default_value(4), "Number of layers")
                ("MESH2D.NodeFile", po::value<std::string>(), "An array of the node coordinates")
                ("MESH2D.Meshfile", po::value<std::string>(), "An array of the Mesh2D ids")
                ("MESH2D.Nfaces", po::value<int>()->default_value(0), "Number of faces per layer")
                ("MESH2D.Nelements", po::value<int>()->default_value(0), "Number of elements per layer")
                ("MESH2D.Nnodes", po::value<int>()->default_value(0), "Number of nodes per layer")
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
            nNodes = vm_vfo["MESH2D.Nnodes"].as<int>();
            nElements = vm_vfo["MESH2D.Nelements"].as<int>();
            nFaces = vm_vfo["MESH2D.Nfaces"].as<int>();
            nTotalFaces = nFaces*nLayers;

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

            if (interp_type == ic::MeshVelInterpType::FACE) {
                std::string faceidsfile = vm_vfo["MESH2D.FaceIdFile"].as<std::string>();
                tf = readFaceIds(faceidsfile);
                if (!tf) { return false;}
                {// read Nodes
                    std::string nodefile = vm_vfo["MESH2D.NodeFile"].as<std::string>();
                    std::vector<std::vector<double>> tmp;
                    tf = READ::read2Darray(nodefile,2,tmp);
                    if (!tf) { return false;}
                    for (unsigned int i = 0; i < tmp.size(); ++i){
                        nds.push_back(vec3(tmp[i][0], tmp[i][1], 0.0));
                    }
                }

                {// read mesh
                    std::string mshfile = vm_vfo["MESH2D.Meshfile"].as<std::string>();
                    std::vector<std::vector<int>> tmp;
                    tf = READ::read2Darray(mshfile,4,tmp);
                    if (!tf) { return false;}
                    for (unsigned int i = 0; i < tmp.size(); ++i){
                        std::vector<int>tmp_ids;
                        for(unsigned int j = 0; j < tmp[i].size(); ++j){
                            if (tmp[i][j] != 0){
                                tmp_ids.push_back(tmp[i][j] - 1);
                            }
                        }
                        msh.push_back(tmp_ids);
                    }
                }
            }
            else if (interp_type == ic::MeshVelInterpType::NODE){
                std::string meshfile = vm_vfo["MESH2D.Meshfile"].as<std::string>();
                tf = readFaceIds(meshfile);
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
                bool tf = readXYZVelocity();
                if (!tf){return false;}
                break;
            }
            case ic::MeshVelInterpType::NODE:
            {
                bool tf = readXYZVelocity();
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
        bool getfid = false;
        std::vector<std::vector<int>> fids;
#if _USEHF > 0
        std::string ext = ic::getExtension(filename);
        if (ext.compare(".h5") == 0){
            const std::string FIDNameSet("FACEID");
            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet datasetFID = HDFNfile.getDataSet(FIDNameSet);
            datasetFID.read(fids);
            getfid = true;
        }
#endif
        if (!getfid){
            bool tf = ic::READ::read2Darray(filename,4,fids);
            if (!tf){return false;}
        }
        for (int i = 0; i < fids.size(); i++){
            std::vector<int> tmp_ids;
            for (int j = 0; j < fids[i].size(); ++j){
                if (fids[i][j] != 0){
                    tmp_ids.push_back(fids[i][j]);
                }
            }
            FaceIds.push_back(tmp_ids);
        }
        return true;
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
                nodeInterpolation(vel, ids, weights, i1, i2, t);
                break;
            }
            case ic::MeshVelInterpType::FACE:
            {
                faceInterpolation(vel, ids, weights, i1, i2, t);
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
        int idx = elid + lay * nElements;
        vel = VEL.getVelocity(idx, i1, i2, t);
    }

    void Mesh2DVel::nodeInterpolation(vec3& vel,
                                      std::vector<int>& ids,
                                      std::vector<double>& weights,
                                      int i1, int i2, double t) {
        int elid = ids[0];
        int lay = ids[1];
        double n1, n2, n3, n4;
        std::vector<double> N;
        if (FaceIds[elid].size() == 4) {
            QuadShapeFunctions(weights[0], weights[1], n1, n2, n3, n4);
            N.push_back(n1);
            N.push_back(n2);
            N.push_back(n3);
            N.push_back(n4);
        }
        else if (FaceIds[elid].size() == 3){
            N.push_back(weights[0]);
            N.push_back(weights[1]);
            N.push_back(1 - weights[0] - weights[1]);
        }


        std::vector<int> idxTop;
        std::vector<int> idxBot;
        for (unsigned int i = 0; i < FaceIds[elid].size(); ++i){
            idxTop.push_back(FaceIds[elid][i]-1 + lay * nNodes);
            idxBot.push_back(FaceIds[elid][i]-1 + (lay+1) * nNodes);
        }
        std::vector<vec3> velTop, velBot;
        VEL.getVelocity(idxTop,i1,i2,t,velTop);
        VEL.getVelocity(idxBot,i1,i2,t,velBot);

        vec3 vt, vb;
        for (unsigned int i = 0; i < velTop.size(); ++i){
            vt = vt + velTop[i]*N[i];
            vb = vb + velBot[i]*N[i];
        }

        vel = vt * weights[2] + vb * (1.0-weights[2]) ;
    }


    void Mesh2DVel::faceInterpolation(ic::vec3& vel,
                                      std::vector<int>& ids,
                                      std::vector<double>& weights,
                                      int i1, int i2, double t) {
        int elid = ids[0];
        int lay = ids[1];

        vec3 tmp;
        int idx;
        if (FaceIds[elid].size() == 4){
            double vf1, vf2, vf3, vf4;// Face velocities
            // face 1
            idx = std::abs(FaceIds[elid][0]) - 1 + lay*nFaces;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf1 = -1.0*sgnFace(FaceIds[elid][0]) * tmp.x;
            // face 2
            idx = std::abs(FaceIds[elid][1]) - 1 + lay*nFaces;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf2 = sgnFace(FaceIds[elid][1]) * tmp.x;
            // face 3
            idx = std::abs(FaceIds[elid][2]) - 1 + lay*nFaces;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf3 = sgnFace(FaceIds[elid][2]) * tmp.x;
            //face 4
            idx = std::abs(FaceIds[elid][3]) - 1 + lay*nFaces;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf4 = -1.0*sgnFace(FaceIds[elid][3]) * tmp.x;

            double u01 = weights[0]/2.0 + 0.5;
            double v01 = weights[1]/2.0 + 0.5;
            double vxt = u01 * vf2 + (1 - u01) * vf4; // Velocity on local space
            double vyt = (1 - v01) * vf1 + v01 * vf3;
            double vlen = std::sqrt(vxt*vxt + vyt*vyt);
            double a, b, c, d;
            calculateQuadJacobian(weights[0], weights[1],
                                  nds[msh[elid][0]], nds[msh[elid][1]],
                                  nds[msh[elid][2]], nds[msh[elid][3]],
                                  a, b, c, d);
            double vx = a * vxt + c * vyt; // Velocity on global space
            double vy = b * vxt + d * vyt;
            double vlen1 = std::sqrt(vx*vx + vy*vy);
            vx = vlen*(vx / vlen1);
            vy = vlen*(vy / vlen1);
            vel.x = vx;
            vel.y = vy;
        }
        else if (FaceIds[elid].size() == 3){
            double vf1, vf2, vf3;// Face velocities
            // face 1
            idx = std::abs(FaceIds[elid][0]) - 1 + lay*nFaces;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf1 = sgnFace(FaceIds[elid][0]) * tmp.x;
            // face 2
            idx = std::abs(FaceIds[elid][1]) - 1 + lay*nFaces;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf2 = sgnFace(FaceIds[elid][1]) * tmp.x;
            // face 3
            idx = std::abs(FaceIds[elid][2]) - 1 + lay*nFaces;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf3 = sgnFace(FaceIds[elid][2]) * tmp.x;

            double u = weights[0];
            double v = weights[1];
            double vxt = u*sqrt2*vf1 + (u-1)*vf2 + u*vf3; // Velocity on local space
            double vyt = v*sqrt2*vf1 + v*vf2 + (v-1)*vf3;
            double vlen = std::sqrt(vxt*vxt + vyt*vyt);
            double a, b, c, d;
            calculateTriangleJacobian(nds[msh[elid][0]], nds[msh[elid][1]],
                                      nds[msh[elid][2]], a, b, c, d);
            double vx = a * vxt + c * vyt; // Velocity on global space
            double vy = b * vxt + d * vyt;
            double vlen1 = std::sqrt(vx*vx + vy*vy);
            // Scale the normalized direction of the global scale with the velocity magnitude
            vx = vlen*(vx / vlen1);
            vy = vlen*(vy / vlen1);
            vel.x = vx;
            vel.y = vy;
        }

        // Top vertical velocity
        idx = nTotalFaces + elid + lay*nElements;
        tmp = VEL.getVelocity(idx, i1, i2, t);
        double vt = tmp.x;
        // Bottom vertical velocity
        idx = nTotalFaces + elid + (lay+1)*nElements;
        tmp = VEL.getVelocity(idx, i1, i2, t);
        double vb = tmp.x;
        vel.z = vt * weights[2] + vb * (1.0-weights[2]);


    }

    bool Mesh2DVel::readXYZVelocity() {
        switch (Vtype) {
            case ic::VelType::STEADY:
            {
                std::string filename = Prefix + ic::num2Padstr(world.rank(), leadingZeros) + Suffix;
                if (Suffix.compare(".h5") == 0){
                    //TODO
                    std::cout << "MESH2D xyz - Steady State - H5 input is not implemented yet" << std::endl;
                    return false;
                }
                std::cout << "\tReading file " + filename << std::endl;
                std::ifstream datafile(filename.c_str());
                if (!datafile.good()){
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
        return false;
    }

    bool Mesh2DVel::readFaceVelocity() {
        bool out = false;
        switch (Vtype){
            case ic::VelType::STEADY:
            {
                std::string filename = Prefix + ic::num2Padstr(world.rank(), leadingZeros) + Suffix;
                if (Suffix.compare(".h5") == 0){
                    //TODO
                    std::cout << "MESH2D for Faces - Steady State - H5 input is not implemented yet" << std::endl;
                    return false;
                }
                std::cout << "\tReading file " + filename << std::endl;
                std::ifstream datafile(filename.c_str());
                if (!datafile.good()){
                    std::cout << "Can't open the file " << filename << std::endl;
                    return false;
                }
                else{
                    std::string line;
                    double vf;
                    std::vector<double> tmp;
                    while (getline(datafile, line)){
                        if (line.size() > 1){
                            std::istringstream inp(line.c_str());
                            inp >> vf;
                            tmp.push_back(vf);
                        }
                    }
                    datafile.close();
                    nPoints = tmp.size();
                    VEL.init(nPoints, nSteps, 1);
                    for (unsigned int i = 0; i < tmp.size(); ++i){
                        VEL.setVELvalue(tmp[i] * multiplier, i, 0, ic::coordDim::vx);
                    }
                    out = true;
                }
                break;
            }
            case ic::VelType::TRANS:
            {
                //TODO
                break;
            }
        }
        return out;
    }

}

#endif //ICHNOS_VELOCITY_MESH2D_H
