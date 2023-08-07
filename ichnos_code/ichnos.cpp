// ichnos.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define _CRT_SECURE_NO_WARNINGS

//const int dbg_rank = 12;

#include <boost/mpi.hpp>

#include <iostream>
#include <ctime>

#include "ichnos_options.h"
#include "ichnos_structures.h"
#include "velocity_base.h"
#include "ichnos_XYZ_base.h"
#include "trace.h"
//#include "iwfmVelocity.h"
//#include "npsatVelocity.h"
#include "stochasticVelocity.h"
#include "TransientVelocity.h"
#include "velocity_mesh2D.h"
#include "XYZ_mesh2D.h"
#include "ichnos_rwpt.h"


ICHNOS::SingletonGenerator* ICHNOS::SingletonGenerator::_instance = nullptr;

int main(int argc, char* argv[])
{
    boost::mpi::environment env( argc, argv );
    boost::mpi::communicator world;

    //{ // Random number generator
        //ICHNOS::SingletonGenerator* RG = RG->getInstance();
        //int cnt = 0;
        //std::vector<int> prop(10,0);
        //for (int i = 0; i < 100000; ++i) {
            //int r = RG->randomNumber(0, prop.size());
         //   float r = RG->randomNumber();
        //    r = RG->randomNumber(0.0f, 1.0f);
        //    int ri = RG->randomNumber(0, 4);
        //    r = RG->randomNumber(0.f, 4.f);
            //prop[r]++;
        //}
        //for (int i = 0; i < prop.size();++i)
        //    std::cout << i << " : " << prop[i] << std::endl;
        //RG->printSeed();
    //}

    if (world.rank() == 0) {
        std::time_t result = std::time(nullptr);
        std::cout << "Ichnos started at " << std::asctime(std::localtime(&result)) << std::endl;
        std::cout << "Ichnos will run using " << world.size() << " processors" << std::endl;
        //std::cout << "Good luck with that!" << std::endl;
        std::cout << "Reading input data..." << std::endl;
    }

    ICHNOS::options OPT(world);
    if (!OPT.readInput(argc, argv))
        return 0;
    if (OPT.Popt.Nthreads > 1 && world.size() > 1){
        std::cout << "You cannot have both: number of threads and number of cores greater than 1 at the same time" << std::endl;
        return 0;
    }

    if (OPT.gatherMode) {
        ICHNOS::gather_particles(OPT);
    }
    else {
        ICHNOS::Domain2D domain(OPT.Dopt);
        bool tf;
        switch (OPT.xyztype){
            case ICHNOS::XYZType::CLOUD:
            {
                world.barrier();
                ICHNOS::XYZ_cloud XYZcloud(world);
                XYZcloud.runAsThread = OPT.Popt.RunAsThread;
                tf = XYZcloud.readXYZdata(OPT.getVelFname());
                world.barrier();
                if (!tf){return 0;}

                switch (OPT.velocityFieldType) {
                    case ICHNOS::VelType::DETRM:
                    {
                        world.barrier();
                        ICHNOS::CloudVel VF(world, XYZcloud);
                        tf = VF.readVelocityField(OPT.getVelFname());
                        world.barrier();
                        if (!tf){return 0;}
                        VF.SetStepOptions(OPT.Popt.StepOpt);

                        ICHNOS::ParticleTrace pt(world, VF, domain, OPT.Popt);
                        world.barrier();
                        if (world.rank() == 0)
                            std::cout << "Tracing particles..." << std::endl;
                        pt.Trace();
                        break;
                    }
                    case ICHNOS::VelType::STOCH:
                    {
                        world.barrier();
                        STOCH::MarkovChainVel MCV(world, XYZcloud);
                        tf = MCV.readVelocityField(OPT.getVelFname());
                        break;
                    }
                    case ICHNOS::VelType::RWPT:
                    {
                        world.barrier();
                        ICHNOS::CloudRWVel VF(world, XYZcloud);
                        tf = VF.readVelocityField(OPT.getVelFname());
                        if (!tf){return 0;}
                        VF.SetStepOptions(OPT.Popt.StepOpt);
                        OPT.Popt.UpdateStepSize = 0;
                        ICHNOS::ParticleTrace pt(world, VF, domain, OPT.Popt);
                        world.barrier();
                        if (world.rank() == 0)
                            std::cout << "Tracing particles..." << std::endl;
                        pt.Trace();
                        break;
                    }
                }
                break;
            }
            case ICHNOS::XYZType::MESH2D:
            {
                ICHNOS::XYZ_MESH2D XYZmesh(world);
                XYZmesh.runAsThread = OPT.Popt.RunAsThread;
                tf = XYZmesh.readXYZdata(OPT.getVelFname());
                world.barrier();
                if (!tf){return 0;}

                switch (OPT.velocityFieldType){
                    case ICHNOS::VelType::DETRM:
                    {
                        world.barrier();
                        ICHNOS::Mesh2DVel VF(world, XYZmesh);
                        tf = VF.readVelocityField(OPT.getVelFname());
                        world.barrier();
                        if (!tf){return 0;}
                        VF.SetStepOptions(OPT.Popt.StepOpt);
                        XYZmesh.SetInterpType(VF.getInterpType());

                        ICHNOS::ParticleTrace pt(world, VF, domain, OPT.Popt);
                        world.barrier();
                        if (world.rank() == 0)
                            std::cout << "Tracing particles..." << std::endl;
                        pt.Trace();
                    }
                    case ICHNOS::VelType::STOCH:
                    {
                        world.barrier();
                        STOCH::MarkovChainVel MCV(world, XYZmesh);
                        tf = MCV.readVelocityField(OPT.getVelFname());

                    }
                }

                break;
            }
       }
    }
    
    //std::cout << world.rank() << " has " << VF.getCloudSize() << " points" << std::endl;
    world.barrier();
    if (world.rank() == 0){
        std::time_t result = std::time(nullptr);
        std::cout << "Ichnos Finished at " << std::asctime(std::localtime(&result)) << std::endl;
    }

    return 0;
   }
