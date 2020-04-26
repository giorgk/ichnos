// ichnos.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <boost/mpi.hpp>

#include <iostream>

#include "ichnos_options.h"
#include "ichnos_structures.h"
#include "velocity_base.h"
#include "trace.h"
#include "iwfmVelocity.h"

ICHNOS::SingletonGenerator* ICHNOS::SingletonGenerator::_instance = nullptr;

int main(int argc, char* argv[])
{
    boost::mpi::environment env{ argc, argv };
    boost::mpi::communicator world;

    //ICHNOS::SingletonGenerator* RG = RG->getInstance();
    //int cnt = 0;
    //std::vector<int> prop(10,0);
    //for (int i = 0; i < 100000; ++i) {
    //    int r = RG->randomNumber(0, prop.size());
    //    prop[r]++;
    //}
    //for (int i = 0; i < prop.size();++i)
    //    std::cout << i << " : " << prop[i] << std::endl;
    //RG->printSeed();

    if (world.rank() == 0) {
        std::cout << "Ichnos will run using " << world.size() << " processors" << std::endl;
        std::cout << "Good luck with that!" << std::endl;

    }

    ICHNOS::options OPT(world);
    if (!OPT.readInput(argc, argv))
        return 0;

    if (OPT.gatherMode) {
        ICHNOS::gather_particles(OPT);
    }
    else {
        switch (OPT.velocityFieldType) {
        case ICHNOS::VelType::Cloud3d:
        {
            ICHNOS::pc3D VF(world);
            VF.readVelocityField(OPT.getVelFname());
            ICHNOS::Domain2D domain(OPT.Dopt);
            ICHNOS::ParticleTrace pt(world, VF, domain, OPT.Popt);
            pt.Trace();
            break;
        }
        case ICHNOS::VelType::IWFM:
        {
            IWFM::iwfmVel VF(world);
            VF.readVelocityField(OPT.getVelFname());
            ICHNOS::Domain2D domain(OPT.Dopt);
            ICHNOS::ParticleTrace pt(world, VF, domain, OPT.Popt);
            for (int i = 0; i < OPT.Popt.Nrealizations; ++i) {
                if (world.rank() == 0)
                    std::cout << "      ====== Realization " << i << "======" << std::endl;
                pt.Trace(i);
            }
            break;
        }
        default:
        {
            std::cout << "Invalid velocity type" << std::endl;
            break;
        }
        }
    }
    
    //std::cout << world.rank() << " has " << VF.getCloudSize() << " points" << std::endl;
    return 0;
   }


//100 1 691636.7 4160195 -100.2
//100 1 646753 4263260 - 50.2
//100 1 854295 3966620 -150.2
