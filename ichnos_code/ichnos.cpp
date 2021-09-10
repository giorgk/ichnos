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
#include "iwfmVelocity.h"
#include "npsatVelocity.h"
#include "stochasticVelocity.h"
#include "TransientVelocity.h"


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
                ICHNOS::XYZ_cloud XYZ(world);
                tf = XYZ.readXYZdata(OPT.getVelFname());
                if (!tf){return 0;}

                switch (OPT.velocityFieldType) {
                    case ICHNOS::VelType::TRANS:
                    {
                        world.barrier();
                        TRANS::transVel VF(world);
                        tf = VF.readVelocityField(OPT.getVelFname(), XYZ.getNpnts());
                        if (!tf){return 0;}
                        VF.SetStepOptions(OPT.Popt.StepOpt);

                        ICHNOS::ParticleTrace pt(world, XYZ, VF, domain, OPT.Popt);
                        world.barrier();
                        if (world.rank() == 0)
                            std::cout << "Tracing particles..." << std::endl;
                        pt.Trace();
                        break;
                    }
                    case ICHNOS::VelType::STOCH:
                    {
                        break;
                    }

                }


                break;
            }
            case ICHNOS::XYZType::IWFM:
            {
                ICHNOS::XYZ_IWFM XYZ(world);
                tf = XYZ.readXYZdata(OPT.getVelFname());
                if (!tf){return 0;}

                switch (OPT.velocityFieldType) {
                    case ICHNOS::VelType::TRANS:
                    {
                        TRANS::transVel VF(world);
                        tf = VF.readVelocityField(OPT.getVelFname(), XYZ.getNpnts());
                        if (!tf){return 0;}
                        VF.SetStepOptions(OPT.Popt.StepOpt);

                        ICHNOS::ParticleTrace pt(world, XYZ, VF, domain, OPT.Popt);
                        if (world.rank() == 0)
                            std::cout << "Tracing particles..." << std::endl;
                        pt.Trace();
                        break;
                    }
                    case ICHNOS::VelType::STOCH:
                    {
                        break;
                    }


                break;
            }
        }

//        switch (OPT.velocityFieldType) {
//        case ICHNOS::VelType::STEADY:
//        {
//            NPSAT::npsatVel VF(world, ICHNOS::VelType::STEADY);
//            VF.readVelocityField(OPT.getVelFname());
//            ICHNOS::Domain2D domain(OPT.Dopt);
//            ICHNOS::ParticleTrace pt(world, VF, domain, OPT.Popt);
//            if (world.rank() == 0)
//                std::cout << "Tracing particles..." << std::endl;
//            pt.Trace();
//            break;
//        }
//        case ICHNOS::VelType::TRANS:
//        {
//            TRANS::transVel VF(world, ICHNOS::VelType::TRANS);
//            VF.readVelocityField(OPT.getVelFname());
//            ICHNOS::Domain2D domain(OPT.Dopt);
//            OPT.Popt.bIsTransient = true;
//            ICHNOS::ParticleTrace pt(world, VF, domain, OPT.Popt);
//            if (world.rank() == 0)
//                std::cout << "Tracing particles..." << std::endl;
//            pt.Trace();
//            break;
//        }
//        case ICHNOS::VelType::STOCH:
//        {
//            STOCH::MarkovChainVel VF(world, ICHNOS::VelType::STOCH);
//            ICHNOS::Domain2D domain(OPT.Dopt);
//            VF.readVelocityField(OPT.getVelFname());
//            ICHNOS::ParticleTrace pt(world, VF, domain, OPT.Popt);
//            for (int i = 0; i < OPT.Popt.Nrealizations; ++i) {
//                if (world.rank() == 0)
//                    std::cout << "||======Realization " << i << "======" << std::endl;
//                pt.Trace(i);
//            }
//            break;
//        }
//        default:
//        {
//            std::cout << "Invalid velocity type" << std::endl;
//            break;
//        }
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


//100 1 691636.7 4160195 -100.2
//100 1 646753 4263260 - 50.2
//100 1 854295 3966620 -150.2
//100 1 804115 4006950 -200.2 Best so far