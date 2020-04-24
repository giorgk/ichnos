// ichnos.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <boost/mpi.hpp>

#include <iostream>

#include "ichnos_options.h"
#include "velocity_base.h"
#include "trace.h"
#include "iwfmVelocity.h"

int main(int argc, char* argv[])
{
    boost::mpi::environment env{ argc, argv };
    boost::mpi::communicator world;

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
        }
        case ICHNOS::VelType::IWFM:
        {
            IWFM::iwfmVel VF(world);
            VF.readVelocityField(OPT.getVelFname());
            ICHNOS::Domain2D domain(OPT.Dopt);
            ICHNOS::ParticleTrace pt(world, VF, domain, OPT.Popt);
            pt.Trace();
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
