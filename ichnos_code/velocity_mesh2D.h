//
// Created by giorg on 11/5/2021.
//

#ifndef ICHNOS_VELOCITY_MESH2D_H
#define ICHNOS_VELOCITY_MESH2D_H

#include "ichnos_structures.h"
#include "velocity_base.h"
#include "ichnos_porosity.h"
#include "ichnos_utils.h"

#if _USEHF > 0
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#endif

namespace ic = ICHNOS;

namespace ICHNOS{
    class Mesh2DVel : public ic::velocityField{
    public:
        Mesh2DVel(boost::mpi::communicator& world_in, XYZ_base &XYZ_in);

        bool readVelocityField(std::string vf_file);
        void calcVelocity(vec3& p, vec3 &vel,
                          std::map<int, double>& proc_map,
                          helpVars& pvlu,
                          bool& out,
                          //std::vector<int> &ids,
                          //std::vector<double> &weights,
                          double time = 0);
        void reset(Streamline& S);
        double stepTimeupdate(helpVars& pvlu);
        void updateStep(helpVars& pvlu);
        void getVec3Data(std::vector<ic::vec3> &data);
        ic::MeshVelInterpType getInterpType(){return interp_type;};

    protected:
        ic::VelTR VEL;
        ic::interpType porType;
        std::vector<std::vector<int>> FaceIds;
        Porosity_base porosity;
        ic::MeshVelInterpType interp_type = ic::MeshVelInterpType::UNKNOWN;
        ic::NodeInterpType VXYnodeInterp = ic::NodeInterpType::TRILINEAR;
        ic::NodeInterpType VZnodeInterp = ic::NodeInterpType::TRILINEAR;

        std::vector<vec3> nds;
        std::vector<std::vector<int>> msh;




        int FrequencyStat;
        int nPoints;
        int nSteps;
        int nLayers;
        // Number of face velocities in the horizontal direction.
        // This is not necessarily the total number of faces.
        int nFaceVelperLayer;
        int nElements;
        int nNodes;
        int nTotalFaces;
        int nTotalHORFaces;
        //ic::TimeData tm_data;

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

    Mesh2DVel::Mesh2DVel(boost::mpi::communicator &world_in, XYZ_base &XYZ_in)
            :
            velocityField(world_in, XYZ_in)
    {
        InterpolateOutsideDomain = false;
    }

    bool Mesh2DVel::readVelocityField(std::string vf_file) {

        if (world.rank() == 0)
            std::cout << "--> Velocity configuration file: " << vf_file << std::endl;

        po::options_description velocityFieldOptions("Velocity field options");
        po::variables_map vm_vfo;
        velocityFieldOptions.add_options()
                // Velocity parameters
                ("Velocity.Prefix", po::value<std::string>(), "Prefix for the filename")
                ("Velocity.LeadingZeros", po::value<int>()->default_value(4), "e.g 0002->4, 000->3")
                ("Velocity.Suffix", po::value<std::string>(), "ending of file after proc id")
                ("Velocity.Type", po::value<std::string>(), "Type of velocity.")
                //("Velocity.IsTransient", po::interpolate<int>()->default_value(0), "0->steady state, 1->Transient state")
                ("Velocity.TimeStepFile", po::value<std::string>(), "This filename with the time steps")
                ("Velocity.TimeInterp", po::value<std::string>(), "Interpolation type between time steps")
                ("Velocity.RepeatTime", po::value<double>()->default_value(0.0), "The number of TimeUnits (e.g. days) to repeat after the and of time steps")
                ("Velocity.Multiplier", po::value<double>()->default_value(1.0), "This is a multiplier to scale velocity")
                //("Velocity.nPoints", po::interpolate<int>()->default_value(0), "This is the number of velocity points")
                ("MESH2D.FaceIdFile", po::value<std::string>(), "Face ids for each element, Required for FACE type")
                ("MESH2D.INTERP", po::value<std::string>(), "Type of interpolation. (ELEMENT, NODE or FACE)")
                ("MESH2D.Nlayers", po::value<int>()->default_value(4), "Number of layers")
                ("MESH2D.NodeFile", po::value<std::string>(), "An array of the node coordinates")
                ("MESH2D.MeshFile", po::value<std::string>(), "An array of the Mesh2D ids")
                ("NODE.VXY", po::value<std::string>(), "TRILINEAR or BILINEAR")
                ("NODE.VZ", po::value<std::string>(), "TRILINEAR or BILINEAR")

                //("MESH2D.Nfaces", po::interpolate<int>()->default_value(0), "Number of faces per layer")
                //("MESH2D.Nelements", po::interpolate<int>()->default_value(0), "Number of elements per layer")
                //("MESH2D.Nnodes", po::interpolate<int>()->default_value(0), "Number of nodes per layer")
                //("Velocity.Scale", po::interpolate<double>()->default_value(1.0), "Scale the domain before velocity calculation")
                //("Velocity.Power", po::interpolate<double>()->default_value(3.0), "Power of the IDW interpolation")
                //("Velocity.InitDiameter", po::interpolate<double>()->default_value(5000), "Initial diameter")
                //("Velocity.InitRatio", po::interpolate<double>()->default_value(1), "Initial ratio")



                // Porosity parameters
                ("Porosity.Value", po::value<std::string>(), "Porosity. Either a file or a single number")

                //General
                ("General.OwnerThreshold", po::value<double>()->default_value(0.75), "Threshold for the processor ownership")
                //("General.Threshold", po::interpolate<double>()->default_value(0.001), "Threshold of distance of IDW")
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
                std::string tp = vm_vfo["MESH2D.INTERP"].as<std::string>();
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
            //nPoints = vm_vfo["Velocity.nPoints"].as<int>();
            nLayers = vm_vfo["MESH2D.Nlayers"].as<int>();
            //nNodes = vm_vfo["MESH2D.Nnodes"].as<int>();
            //nElements = vm_vfo["MESH2D.Nelements"].as<int>();
            //nFaceVels = vm_vfo["MESH2D.Nfaces"].as<int>();
            //nTotalFaces = nFaceVels*nLayers;

            //isVeltrans  = vm_vfo["Velocity.IsTransient"].as<int>() != 0;

            // Read the time step
            std::vector<double> TimeSteps;
            isVeltrans = false;
            if (vm_vfo.count("Velocity.TimeStepFile")){
                std::string TSfile = vm_vfo["Velocity.TimeStepFile"].as<std::string>();
                if (!TSfile.empty()){
                    isVeltrans = true;
                }
            }

            if (isVeltrans){
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

            nSteps = TimeSteps.size()-1;

            Prefix = vm_vfo["Velocity.Prefix"].as<std::string>();
            std::string suffix;
            if (vm_vfo.count("Velocity.Suffix")) {
                Suffix = vm_vfo["Velocity.Suffix"].as<std::string>();
            }
            else {
                Suffix = ".ich";
            }
            leadingZeros = vm_vfo["Velocity.LeadingZeros"].as<int>();

            bool tf = readVelocityFiles();

            if (world.size() > 1 && XYZ.runAsThread == 0){
                std::cout << "Proc: " << world.rank() << " has : " << VEL.getNpoints() << "  Velocities "  << std::endl;
            }

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
                    std::string mshfile = vm_vfo["MESH2D.MeshFile"].as<std::string>();
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
                std::string meshfile = vm_vfo["MESH2D.MeshFile"].as<std::string>();
                tf = readFaceIds(meshfile);

                if (vm_vfo.count("NODE.VXY")){
                    std::string vxyinterp = vm_vfo["NODE.VXY"].as<std::string>();
                    if (vxyinterp.compare("BILINEAR") == 0){
                        VXYnodeInterp = ic::NodeInterpType::BILINEAR;
                    }
                }

                if (vm_vfo.count("NODE.VZ")){
                    std::string vzinterp = vm_vfo["NODE.VZ"].as<std::string>();
                    if (vzinterp.compare("BILINEAR") == 0){
                        VZnodeInterp = ic::NodeInterpType::BILINEAR;
                    }
                }

                if (!tf) { return false;}
            }
        }

        { // Porosity
            if (vm_vfo.count("Porosity.Value")){
                std::string porfile = vm_vfo["Porosity.Value"].as<std::string>();
                bool tf = porosity.readData(porfile);
                if (!tf)
                    return false;
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
                nElements = XYZ.getINTInfo(infoType::Nelem);
                break;
            }
            case ic::MeshVelInterpType::NODE:
            {
                bool tf = readXYZVelocity();
                if (!tf){return false;}
                nElements = XYZ.getINTInfo(infoType::Nelem);
                nNodes = XYZ.getINTInfo(infoType::NNodes);
                break;
            }
            case ic::MeshVelInterpType::FACE:
            {
                nElements = XYZ.getINTInfo(infoType::Nelem);
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

    void Mesh2DVel::calcVelocity(vec3& p, vec3 &vel,
                                 std::map<int, double>& proc_map,
                                 helpVars& pvlu,
                                 bool& out,
                                 //std::vector<int> &ids,
                                 //std::vector<double> &weights,
                                 double tm) {
        std::vector<int> ids;
        std::vector<double> weights;
        XYZ.calcWeights(p, ids,weights, proc_map, pvlu, out);
        if (!out){
            vel.makeInvalid();
            return;
        }


        if (ids[0] == -9 || ids[1] == -9){
            vel.makeInvalid();
            return;
        }
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

        double por = porosity.interpolate(p);
        vel = vel * (1/por);

        //tm_data.tm = VEL.getTSvalue(i1) * (1-t) + VEL.getTSvalue(i2)*t;
        //tm_data.idx1 = i1;
        //tm_data.idx2 = i2;
        //tm_data.t = t;
    }

    double Mesh2DVel::stepTimeupdate(ICHNOS::helpVars &pvlu) {
        return VEL.stepTimeUpdate(pvlu, stepOpt);
    }

    void Mesh2DVel::updateStep(helpVars& pvlu) {

    }

    void Mesh2DVel::getVec3Data(std::vector<ic::vec3>& data) {

    }

    void Mesh2DVel::reset(Streamline& S) {
        porosity.reset();
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
            if (VXYnodeInterp == ic::NodeInterpType::BILINEAR){
                // We will replace the x y velocities of the top layer to the bottom layer
                velBot[i].x = velTop[i].x;
                velBot[i].y = velTop[i].y;
            }
            if (VZnodeInterp == ic::NodeInterpType::BILINEAR){
                velBot[i].z = velTop[i].z;
            }

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
            idx = std::abs(FaceIds[elid][0]) - 1 + lay * nFaceVelperLayer;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf1 = -1.0*sgnFace(FaceIds[elid][0]) * tmp.x;
            // face 2
            idx = std::abs(FaceIds[elid][1]) - 1 + lay * nFaceVelperLayer;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf2 = sgnFace(FaceIds[elid][1]) * tmp.x;
            // face 3
            idx = std::abs(FaceIds[elid][2]) - 1 + lay * nFaceVelperLayer;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf3 = sgnFace(FaceIds[elid][2]) * tmp.x;
            //face 4
            idx = std::abs(FaceIds[elid][3]) - 1 + lay * nFaceVelperLayer;
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
            idx = std::abs(FaceIds[elid][0]) - 1 + lay * nFaceVelperLayer;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf1 = sgnFace(FaceIds[elid][0]) * tmp.x;
            // face 2
            idx = std::abs(FaceIds[elid][1]) - 1 + lay * nFaceVelperLayer;
            tmp = VEL.getVelocity(idx, i1, i2, t);
            vf2 = sgnFace(FaceIds[elid][1]) * tmp.x;
            // face 3
            idx = std::abs(FaceIds[elid][2]) - 1 + lay * nFaceVelperLayer;
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
        idx = nTotalHORFaces + elid + lay*nElements;
        tmp = VEL.getVelocity(idx, i1, i2, t);
        double vt = tmp.x;
        // Bottom vertical velocity
        idx = nTotalHORFaces + elid + (lay+1)*nElements;
        tmp = VEL.getVelocity(idx, i1, i2, t);
        double vb = tmp.x;
        vel.z = vt * weights[2] + vb * (1.0-weights[2]);


    }

    bool Mesh2DVel::readXYZVelocity() {

        int proc_id = world.rank();
        if (XYZ.runAsThread){
            proc_id = 0;
        }

        if (!isVeltrans){
            std::string filename = Prefix + ic::num2Padstr(proc_id, leadingZeros) + Suffix;

            if (Suffix.compare(".h5") == 0){
#if _USEHF > 0
                std::vector<std::vector<double>> VXYZ;
                bool tf = READ::H5SteadyState3DVelocity(filename, VXYZ);
                if (tf){
                    nPoints = VXYZ[0].size();
                    VEL.init(nPoints,1, 3);
                    for (int i = 0; i < VXYZ[0].size(); i++){
                        VEL.setVELvalue(VXYZ[0][i]*multiplier, i, 0, coordDim::vx);
                        VEL.setVELvalue(VXYZ[1][i]*multiplier, i, 0, coordDim::vy);
                        VEL.setVELvalue(VXYZ[2][i]*multiplier, i, 0, coordDim::vz);
                    }
                    return true;
                }
#endif
            }
            else{
                std::vector<std::vector<double>> data;
                bool tf = READ::read2Darray<double>(filename, 3, data);
                if (tf){
                    nPoints = static_cast<int>(data.size());
                    VEL.init(nPoints, 1, 3);
                    for (int i = 0; i < nPoints; i++){
                        VEL.setVELvalue(data[i][0]*multiplier, i, 0, coordDim::vx);
                        VEL.setVELvalue(data[i][1]*multiplier, i, 0, coordDim::vy);
                        VEL.setVELvalue(data[i][2]*multiplier, i, 0, coordDim::vz);
                    }
                    return true;
                }
            }
        }
        else{
            if (Suffix.compare(".h5") == 0){
#if _USEHF > 0
                std::string filename = Prefix + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
                bool tf = READ::H5Transient3DVelocity(filename, nSteps, multiplier, VEL);
                return tf;
#endif
            }
            else{
                std::string filenameVX = Prefix + "VX_" + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
                std::string filenameVY = Prefix + "VY_" + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
                std::string filenameVZ = Prefix + "VZ_" + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
                std::vector<std::vector<double>> VX, VY, VZ;
                bool tf = READ::read2Darray<double>(filenameVX, nSteps, VX);
                if (!tf)
                    return false;
                tf = READ::read2Darray<double>(filenameVY, nSteps, VY);
                if (!tf)
                    return false;
                tf = READ::read2Darray<double>(filenameVZ, nSteps, VZ);
                if (!tf)
                    return false;

                if (VX[0].size() != nSteps){
                    std::cout << "The number of time steps (" << VX[0].size() << ") in the file " << filenameVX
                              << "\n does not match the number of steps in the Time step file " << nSteps << std::endl;
                    return false;
                }
                VEL.init(VX.size(), nSteps, 3);
                for (int i = 0; i < VX.size(); i++){
                    for (int j = 0; j < VX[i].size(); ++j) {
                        VEL.setVELvalue(VX[i][j]*multiplier, i, j, ICHNOS::coordDim::vx);
                        VEL.setVELvalue(VY[i][j]*multiplier, i, j, ICHNOS::coordDim::vy);
                        VEL.setVELvalue(VZ[i][j]*multiplier, i, j, ICHNOS::coordDim::vz);
                    }
                }
                return true;
            }
        }
        return false;
    }

    bool Mesh2DVel::readFaceVelocity() {
        bool out = false;
        bool bisascii = true;
        int proc_id = world.rank();
        if (XYZ.runAsThread){
            proc_id = 0;
        }

        auto start = std::chrono::high_resolution_clock::now();

        if (!isVeltrans){
            std::vector<std::vector<double>> VHOR;
            std::vector<std::vector<double>> VVER;
            if (Suffix.compare(".h5") == 0){
#if _USEHF > 0
                //TODO
                std::string filename = Prefix + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
                const std::string VHORNameSet("VHOR");
                const std::string VVERNameSet("VVER");
                HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
                HighFive::DataSet datasetVHOR = HDFNfile.getDataSet(VHORNameSet);
                HighFive::DataSet datasetVVER = HDFNfile.getDataSet(VVERNameSet);

                datasetVHOR.read(VHOR);
                datasetVVER.read(VVER);
                bisascii = false;
#endif
            }
            else{
                std::string filenameHOR = Prefix + "VHOR_" + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
                std::string filenameVER = Prefix + "VVER_" + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
                bool tf = READ::read2Darray<double>(filenameHOR, 1, VHOR);
                if (!tf){
                    return false;
                }
                tf = READ::read2Darray<double>(filenameVER, 1, VVER);
                if (!tf){
                    return false;
                }
            }

            int nTotalVERface;
            if (bisascii){
                nTotalHORFaces = static_cast<int>(VHOR.size());
                nTotalVERface = static_cast<int>(VVER.size());
            }
            else{
                nTotalHORFaces = static_cast<int>(VHOR[0].size());
                nTotalVERface = static_cast<int>(VVER[0].size());
            }

            nFaceVelperLayer = nTotalHORFaces/nLayers;
            nTotalFaces = nTotalHORFaces + nTotalVERface;
            int n = nTotalVERface / (nLayers + 1);
            if (nElements != n){
                std::cout << "The number of elements in Element file is not the same as the elements in the VVER file" << std::endl;
                return false;
            }
            VEL.init(nTotalFaces, 1, 1);

            for (int i = 0; i < nTotalHORFaces; ++i){
                if (bisascii){
                    VEL.setVELvalue(VHOR[i][0] * multiplier, i, 0, ic::coordDim::vx);
                }
                else{
                    VEL.setVELvalue(VHOR[0][i] * multiplier, i, 0, ic::coordDim::vx);
                }
            }

            for (int i = 0; i < nTotalVERface; ++i){
                if (bisascii){
                    VEL.setVELvalue(VVER[i][0] * multiplier, i + nTotalHORFaces, 0, ic::coordDim::vx);
                }
                else{
                    VEL.setVELvalue(VVER[0][i] * multiplier, i + nTotalHORFaces, 0, ic::coordDim::vx);
                }
            }
            out = true;

        }// If velocity is transient
        else{
            std::vector<std::vector<double>> VHOR;
            std::vector<std::vector<double>> VVER;
            if (Suffix.compare(".h5") == 0) {
#if _USEHF > 0
                // TODO
                std::string filename = Prefix + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
                const std::string VHORNameSet("VHOR");
                const std::string VVERNameSet("VVER");
                HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
                HighFive::DataSet datasetVHOR = HDFNfile.getDataSet(VHORNameSet);
                HighFive::DataSet datasetVVER = HDFNfile.getDataSet(VVERNameSet);

                datasetVHOR.read(VHOR);
                datasetVVER.read(VVER);
                bisascii = false;
#endif
            }
            else{
                std::string filenameHOR = Prefix + "VHOR_" + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
                std::string filenameVER = Prefix + "VVER_" + ic::num2Padstr(proc_id, leadingZeros) + Suffix;
                bool tf = READ::read2Darray<double>(filenameHOR, nSteps, VHOR);
                if (VHOR[0].size() != nSteps){
                    std::cout << "The number of time steps (" << VHOR[0].size() << ") in the file " << filenameHOR
                              << "\n does not match the number of steps in the Time step file " << nSteps << std::endl;
                    return false;
                }
                if (!tf){
                    return false;
                }
                tf = READ::read2Darray<double>(filenameVER, nSteps, VVER);
                if (!tf){
                    return false;
                }
            }

            int nTotalVERface;
            if (bisascii){
                nTotalHORFaces = static_cast<int>(VHOR.size());
                nTotalVERface = static_cast<int>(VVER.size());
            }
            else{
                nTotalHORFaces = static_cast<int>(VHOR[0].size());
                nTotalVERface = static_cast<int>(VVER[0].size());
            }

            nFaceVelperLayer = nTotalHORFaces/nLayers;
            nTotalFaces = nTotalHORFaces + nTotalVERface;
            int n = nTotalVERface / (nLayers + 1);
            if (nElements != n){
                bisascii = false;
                std::cout << "The number of elements in Element file is not the same as the elements in the VVER file" << std::endl;
                return false;
            }
            VEL.init(nTotalFaces, nSteps, 1);
            for (int i = 0; i < nTotalHORFaces; ++i){
                for (int j = 0; j < nSteps; ++j){
                    if (bisascii){
                        VEL.setVELvalue(VHOR[i][j] * multiplier, i, j, ic::coordDim::vx);
                    }
                    else{
                        VEL.setVELvalue(VHOR[j][i] * multiplier, i, j, ic::coordDim::vx);
                    }
                }
            }

            for (int i = 0; i < nTotalVERface; ++i){
                for (int j = 0; j < nSteps; ++j){
                    if (bisascii){
                        VEL.setVELvalue(VVER[i][j] * multiplier, i + nTotalHORFaces, j, ic::coordDim::vx);
                    }
                    else{
                        VEL.setVELvalue(VVER[j][i] * multiplier, i + nTotalHORFaces, j, ic::coordDim::vx);
                    }
                }
            }
            out = true;
        }

        if (out){
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            std::cout << "\tRead Velocity in : " << elapsed.count() << std::endl;
        }

        return out;
    }

}

#endif //ICHNOS_VELOCITY_MESH2D_H
